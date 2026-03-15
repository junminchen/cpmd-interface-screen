#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import json
import re
import struct
from dataclasses import dataclass
from pathlib import Path

import numpy as np
import openmm.app as app


E_CHARGE_C = 1.602176634e-19
ANGSTROM2_TO_M2 = 1.0e-20
F_PER_M2_TO_UF_PER_CM2 = 100.0
DEFAULT_DONOR_ELEMENTS = {"O", "N", "S"}


def _read_record(handle):
    head = handle.read(4)
    if len(head) == 0:
        return None
    if len(head) < 4:
        raise EOFError("Truncated DCD record head")
    n = struct.unpack("<i", head)[0]
    payload = handle.read(n)
    if len(payload) != n:
        raise EOFError("Truncated DCD record payload")
    tail = handle.read(4)
    if len(tail) < 4:
        raise EOFError("Truncated DCD record tail")
    if struct.unpack("<i", tail)[0] != n:
        raise ValueError("DCD record length mismatch")
    return payload


def dcd_frame_iterator(dcd_path: Path):
    with dcd_path.open("rb") as handle:
        rec = _read_record(handle)
        if rec is None:
            return
        if len(rec) < 4 or rec[:4] != b"CORD":
            raise ValueError("Unsupported DCD header (expected CORD)")

        rec = _read_record(handle)
        if rec is None:
            raise ValueError("Missing DCD title block")

        rec = _read_record(handle)
        if rec is None or len(rec) != 4:
            raise ValueError("Missing DCD natoms block")
        natom = struct.unpack("<i", rec)[0]

        while True:
            rec = _read_record(handle)
            if rec is None:
                break
            if len(rec) == 48:
                rec = _read_record(handle)
                if rec is None:
                    break

            if len(rec) != 4 * natom:
                raise ValueError("Unexpected DCD x block length")
            x = np.frombuffer(rec, dtype=np.float32, count=natom)

            recy = _read_record(handle)
            recz = _read_record(handle)
            if recy is None or recz is None:
                raise EOFError("Truncated DCD coordinate blocks")
            if len(recy) != 4 * natom or len(recz) != 4 * natom:
                raise ValueError("Unexpected DCD y/z block length")

            y = np.frombuffer(recy, dtype=np.float32, count=natom)
            z = np.frombuffer(recz, dtype=np.float32, count=natom)
            yield x, y, z


def resolve_input_path(raw: str, script_dir: Path) -> Path:
    candidate = Path(raw).expanduser()
    if candidate.is_absolute():
        return candidate
    cwd_path = (Path.cwd() / candidate).resolve()
    if cwd_path.exists():
        return cwd_path
    return (script_dir / candidate).resolve()


def resolve_output_path(raw: str, script_dir: Path) -> Path:
    candidate = Path(raw).expanduser()
    if candidate.is_absolute():
        return candidate
    return (Path.cwd() / candidate if Path.cwd().joinpath(candidate).parent.exists() else script_dir / candidate).resolve()


def slugify(text: str) -> str:
    text = re.sub(r"[^0-9a-zA-Z]+", "_", text.strip().lower()).strip("_")
    return text or "species"


def parse_charge_log(path: Path):
    steps = []
    q_cath = []
    q_anode = []
    q_total = []
    for line in path.read_text().splitlines():
        line = line.strip()
        if (not line) or line.startswith("#"):
            continue
        toks = line.split()
        if len(toks) < 4:
            continue
        steps.append(int(float(toks[0])))
        q_cath.append(float(toks[1]))
        q_anode.append(float(toks[2]))
        q_total.append(float(toks[3]))
    if not steps:
        raise RuntimeError(f"No data rows in {path}")
    return np.asarray(steps), np.asarray(q_cath), np.asarray(q_anode), np.asarray(q_total)


def block_stats(values: np.ndarray, nblocks: int):
    n = values.shape[0]
    if n < 2:
        return float(values.mean()), float("nan"), 1
    bsize = max(1, n // max(1, nblocks))
    nblk = max(1, n // bsize)
    if nblk < 2:
        return float(values.mean()), float("nan"), nblk
    means = []
    for block in range(nblk):
        lo = block * bsize
        hi = min((block + 1) * bsize, n)
        if hi <= lo:
            continue
        means.append(float(np.mean(values[lo:hi])))
    if len(means) < 2:
        return float(values.mean()), float("nan"), len(means)
    arr = np.asarray(means)
    return float(arr.mean()), float(arr.std(ddof=1)), len(means)


@dataclass
class SpeciesDefinition:
    name: str
    slug: str
    resnames: set[str]
    donor_atoms: set[str] | None
    donor_elements: set[str]


@dataclass
class SpeciesSelection:
    definition: SpeciesDefinition
    residue_atom_indices: list[np.ndarray]
    residue_masses: list[np.ndarray]
    donor_atom_indices: np.ndarray


def load_species_definitions(items: list[dict], label: str) -> list[SpeciesDefinition]:
    definitions = []
    for item in items:
        name = str(item["name"]).strip()
        if not name:
            raise ValueError(f"{label}.name cannot be empty")
        resnames = {str(v).strip() for v in item["resnames"]}
        donor_atoms = item.get("donor_atoms")
        donor_atoms_set = {str(v).strip() for v in donor_atoms} if donor_atoms else None
        donor_elements = item.get("donor_elements")
        donor_elements_set = (
            {str(v).strip().upper() for v in donor_elements}
            if donor_elements
            else set(DEFAULT_DONOR_ELEMENTS)
        )
        definitions.append(
            SpeciesDefinition(
                name=name,
                slug=slugify(name),
                resnames=resnames,
                donor_atoms=donor_atoms_set,
                donor_elements=donor_elements_set,
            )
        )
    return definitions


def atom_matches_donor(atom, definition: SpeciesDefinition) -> bool:
    if definition.donor_atoms is not None:
        return atom.name in definition.donor_atoms
    if atom.element is None:
        return False
    return atom.element.symbol.upper() in definition.donor_elements


def select_species(topology: app.Topology, masses: np.ndarray, definitions: list[SpeciesDefinition], strict: bool = True):
    selections: dict[str, SpeciesSelection] = {}
    for definition in definitions:
        residue_atom_indices = []
        residue_masses = []
        donor_atom_indices = []
        matched_residues = 0
        for residue in topology.residues():
            if residue.name.strip() not in definition.resnames:
                continue
            matched_residues += 1
            atom_indices = np.asarray([atom.index for atom in residue.atoms()], dtype=int)
            residue_atom_indices.append(atom_indices)
            residue_masses.append(masses[atom_indices])
            donor_atom_indices.extend(
                atom.index for atom in residue.atoms() if atom_matches_donor(atom, definition)
            )
        if matched_residues == 0 and strict:
            names = ", ".join(sorted(definition.resnames))
            raise RuntimeError(f"Did not find residues for species '{definition.name}' ({names})")
        donor_array = np.asarray(donor_atom_indices, dtype=int)
        if matched_residues > 0 and donor_array.size == 0 and strict:
            raise RuntimeError(f"Species '{definition.name}' matched residues but no donor atoms were selected")
        selections[definition.slug] = SpeciesSelection(
            definition=definition,
            residue_atom_indices=residue_atom_indices,
            residue_masses=residue_masses,
            donor_atom_indices=donor_array,
        )
    return selections


def find_lithium_indices(topology: app.Topology, lithium_resnames: set[str]) -> np.ndarray:
    indices = []
    for residue in topology.residues():
        if residue.name.strip() not in lithium_resnames:
            continue
        for atom in residue.atoms():
            if atom.element is not None and atom.element.symbol.upper() == "LI":
                indices.append(atom.index)
        if not indices:
            indices.extend(atom.index for atom in residue.atoms())
    if not indices:
        names = ", ".join(sorted(lithium_resnames))
        raise RuntimeError(f"Did not find Li atoms for residues: {names}")
    return np.asarray(indices, dtype=int)


def compute_residue_com_z(z_coords: np.ndarray, residue_atom_indices: list[np.ndarray], residue_masses: list[np.ndarray]) -> np.ndarray:
    out = np.zeros(len(residue_atom_indices), dtype=float)
    for idx, (atom_idx, mass) in enumerate(zip(residue_atom_indices, residue_masses)):
        out[idx] = float(np.dot(z_coords[atom_idx], mass) / np.sum(mass))
    return out


def donor_contact_counts(
    li_coords: np.ndarray,
    donor_coords: np.ndarray,
    box_x_angstrom: float,
    box_y_angstrom: float,
    use_xy_pbc: bool,
    cutoff_sq: float,
) -> np.ndarray:
    if li_coords.size == 0 or donor_coords.size == 0:
        return np.zeros(li_coords.shape[0], dtype=float)
    dx = donor_coords[None, :, 0] - li_coords[:, None, 0]
    dy = donor_coords[None, :, 1] - li_coords[:, None, 1]
    dz = donor_coords[None, :, 2] - li_coords[:, None, 2]
    if use_xy_pbc:
        dx -= box_x_angstrom * np.round(dx / box_x_angstrom)
        dy -= box_y_angstrom * np.round(dy / box_y_angstrom)
    dist_sq = dx * dx + dy * dy + dz * dz
    return np.sum(dist_sq <= cutoff_sq, axis=1, dtype=float)


def density_from_counts(total_count: float, nframes: int, area_nm2: float, width_angstrom: float) -> float:
    width_nm = width_angstrom * 0.1
    if nframes <= 0 or area_nm2 <= 0.0 or width_nm <= 0.0:
        return float("nan")
    return float(total_count / (nframes * area_nm2 * width_nm))


def summarize_events(event_lengths_frames: list[int], dt_ps: float) -> tuple[float, float, int]:
    if not event_lengths_frames:
        return 0.0, 0.0, 0
    arr = np.asarray(event_lengths_frames, dtype=float)
    return float(arr.mean() * dt_ps), float(np.percentile(arr, 90.0) * dt_ps), int(arr.size)


def safe_ratio(num: float, den: float) -> float:
    if den == 0.0 or not np.isfinite(den):
        return float("nan")
    return float(num / den)


def maybe_compute_capacitance(config: dict, charge_log_path: Path | None, skip_rows: int, nblocks: int):
    if charge_log_path is None or not charge_log_path.exists():
        return {}
    steps, q_cath, q_anode, _ = parse_charge_log(charge_log_path)
    if skip_rows > 0:
        if skip_rows >= len(steps):
            raise RuntimeError(f"--capacitance-skip {skip_rows} >= samples {len(steps)}")
        q_cath = q_cath[skip_rows:]
        q_anode = q_anode[skip_rows:]
    area_A2 = float(config["cell"]["box_x_angstrom"]) * float(config["cell"]["box_y_angstrom"])
    area_m2 = area_A2 * ANGSTROM2_TO_M2
    voltage_v = float(config["electrode"]["voltage_v"])
    c_cath = np.abs(q_cath) * E_CHARGE_C / (area_m2 * voltage_v)
    c_anode = np.abs(q_anode) * E_CHARGE_C / (area_m2 * voltage_v)
    c_avg = 0.5 * (c_cath + c_anode)
    mean_avg, std_avg, blocks = block_stats(c_avg, nblocks)
    return {
        "interfacial_capacitance_uF_per_cm2": mean_avg * F_PER_M2_TO_UF_PER_CM2,
        "interfacial_capacitance_block_std_uF_per_cm2": std_avg * F_PER_M2_TO_UF_PER_CM2
        if np.isfinite(std_avg)
        else float("nan"),
        "interfacial_capacitance_blocks_used": blocks,
    }


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Convert constant-potential MD outputs into additive-screening descriptors."
    )
    parser.add_argument("--top", required=True, help="Topology PDB")
    parser.add_argument("--traj", required=True, help="Trajectory DCD")
    parser.add_argument("--config", required=True, help="MD config JSON with cell/md settings")
    parser.add_argument("--descriptor-config", required=True, help="Descriptor definition JSON")
    parser.add_argument("--charge-log", default=None, help="Optional electrode charge log for capacitance")
    parser.add_argument("--out-csv", required=True, help="Output CSV with one descriptor row")
    parser.add_argument("--out-json", required=True, help="Output JSON with full descriptor payload")
    parser.add_argument("--formulation-id", default=None, help="Stable formulation identifier for later label joins")
    parser.add_argument("--stride", type=int, default=1, help="Analyze every Nth trajectory frame")
    parser.add_argument("--max-frames", type=int, default=None, help="Stop after this many analyzed frames")
    parser.add_argument("--frame-dt-ps", type=float, default=None, help="Override time spacing between stored DCD frames")
    parser.add_argument("--capacitance-skip", type=int, default=0, help="Rows to skip from charge log")
    parser.add_argument("--capacitance-nblocks", type=int, default=5)
    args = parser.parse_args()

    script_dir = Path(__file__).resolve().parent
    top_path = resolve_input_path(args.top, script_dir)
    traj_path = resolve_input_path(args.traj, script_dir)
    config_path = resolve_input_path(args.config, script_dir)
    descriptor_config_path = resolve_input_path(args.descriptor_config, script_dir)
    charge_log_path = resolve_input_path(args.charge_log, script_dir) if args.charge_log else None
    out_csv_path = resolve_output_path(args.out_csv, script_dir)
    out_json_path = resolve_output_path(args.out_json, script_dir)

    config = json.loads(config_path.read_text())
    descriptor_config = json.loads(descriptor_config_path.read_text())

    pdb = app.PDBFile(str(top_path))
    topology = pdb.topology
    atoms = list(topology.atoms())
    natom = len(atoms)
    masses = np.asarray(
        [(atom.element.mass._value if atom.element is not None else 0.0) for atom in atoms],
        dtype=float,
    )
    pdb_positions_angstrom = np.asarray([[v.x, v.y, v.z] for v in pdb.positions], dtype=float) * 10.0

    lithium_resnames = {str(v).strip() for v in descriptor_config["lithium_resnames"]}
    additive_defs = load_species_definitions(descriptor_config["additive_species"], "additive_species")
    reference_defs = load_species_definitions(descriptor_config["reference_species"], "reference_species")

    additive_sel = select_species(topology, masses, additive_defs, strict=True)
    reference_sel = select_species(topology, masses, reference_defs, strict=True)
    li_atom_indices = find_lithium_indices(topology, lithium_resnames)

    box_x_angstrom = float(config["cell"]["box_x_angstrom"])
    box_y_angstrom = float(config["cell"]["box_y_angstrom"])
    area_nm2 = (box_x_angstrom * 0.1) * (box_y_angstrom * 0.1)

    z_cathode = descriptor_config.get("z_cathode_surface_angstrom", config["cell"].get("z_cathode_surface_angstrom"))
    z_anode = descriptor_config.get("z_anode_surface_angstrom", config["cell"].get("z_anode_surface_angstrom"))
    if z_cathode is None or z_anode is None:
        chains = list(topology.chains())
        if len(chains) < 2:
            raise RuntimeError("Could not infer electrode surfaces: need config cell surfaces or two electrode chains")
        cath_idx = np.asarray([atom.index for atom in chains[0].atoms()], dtype=int)
        anode_idx = np.asarray([atom.index for atom in chains[1].atoms()], dtype=int)
        z_cathode = float(np.max(pdb_positions_angstrom[cath_idx, 2]))
        z_anode = float(np.min(pdb_positions_angstrom[anode_idx, 2]))
    else:
        z_cathode = float(z_cathode)
        z_anode = float(z_anode)

    electrode_margin = float(descriptor_config.get("electrode_margin_angstrom", 1.0))
    interface_width = float(descriptor_config.get("interface_width_angstrom", 5.0))
    cath_lo = z_cathode + electrode_margin
    cath_hi = cath_lo + interface_width
    anode_hi = z_anode - electrode_margin
    anode_lo = anode_hi - interface_width

    bulk_cfg = descriptor_config.get("bulk_region", {})
    if bulk_cfg.get("mode", "auto") == "auto":
        bulk_zmin = cath_hi
        bulk_zmax = anode_lo
    else:
        bulk_zmin = float(bulk_cfg["zmin_angstrom"])
        bulk_zmax = float(bulk_cfg["zmax_angstrom"])
    if not (bulk_zmax > bulk_zmin):
        raise RuntimeError(
            f"Invalid bulk region [{bulk_zmin:.3f}, {bulk_zmax:.3f}); reduce interface width or override bulk_region"
        )

    report_interval = float(config["md"]["report_interval"])
    timestep_fs = float(config["md"]["timestep_fs"])
    frame_dt_ps = float(args.frame_dt_ps if args.frame_dt_ps is not None else timestep_fs * report_interval * 1.0e-3)
    effective_dt_ps = frame_dt_ps * max(1, args.stride)
    li_cutoff = float(descriptor_config.get("li_coordination_cutoff_angstrom", 2.8))
    cutoff_sq = li_cutoff * li_cutoff
    use_xy_pbc = bool(descriptor_config.get("use_xy_pbc", True))

    additive_stats = {}
    for slug, selection in additive_sel.items():
        nmol = len(selection.residue_atom_indices)
        additive_stats[slug] = {
            "name": selection.definition.name,
            "n_molecules": nmol,
            "cathode_count_sum": 0.0,
            "anode_count_sum": 0.0,
            "bulk_count_sum": 0.0,
            "occupancy_sum": 0.0,
            "molecule_frame_count": 0.0,
            "current_run_lengths": np.zeros(nmol, dtype=int),
            "event_lengths_frames": [],
            "li_cn_sum": {"all": 0.0, "interface": 0.0, "bulk": 0.0},
            "li_contact_count": {"all": 0.0, "interface": 0.0, "bulk": 0.0},
        }

    reference_cn_sum = {"all": 0.0, "interface": 0.0, "bulk": 0.0}
    li_samples = {"all": 0.0, "interface": 0.0, "bulk": 0.0}
    analyzed_frames = 0

    for iframe, (x, y, z) in enumerate(dcd_frame_iterator(traj_path), start=1):
        if args.stride > 1 and ((iframe - 1) % args.stride != 0):
            continue
        if args.max_frames is not None and analyzed_frames >= args.max_frames:
            break
        if z.shape[0] != natom:
            raise ValueError(f"DCD natom mismatch: frame {iframe} has {z.shape[0]}, topology has {natom}")

        analyzed_frames += 1

        coords = np.column_stack((x.astype(float), y.astype(float), z.astype(float)))

        for slug, selection in additive_sel.items():
            stats = additive_stats[slug]
            z_com = compute_residue_com_z(z, selection.residue_atom_indices, selection.residue_masses)
            in_cathode = (z_com >= cath_lo) & (z_com < cath_hi)
            in_anode = (z_com >= anode_lo) & (z_com < anode_hi)
            in_bulk = (z_com >= bulk_zmin) & (z_com < bulk_zmax)
            in_first_layer = in_cathode | in_anode

            stats["cathode_count_sum"] += float(np.count_nonzero(in_cathode))
            stats["anode_count_sum"] += float(np.count_nonzero(in_anode))
            stats["bulk_count_sum"] += float(np.count_nonzero(in_bulk))
            stats["occupancy_sum"] += float(np.count_nonzero(in_first_layer))
            stats["molecule_frame_count"] += float(z_com.size)

            runs = stats["current_run_lengths"]
            runs[in_first_layer] += 1
            ended = (~in_first_layer) & (runs > 0)
            if np.any(ended):
                stats["event_lengths_frames"].extend(runs[ended].tolist())
                runs[ended] = 0

        li_coords = coords[li_atom_indices]
        li_z = li_coords[:, 2]
        li_region_masks = {
            "all": np.ones(li_coords.shape[0], dtype=bool),
            "interface": ((li_z >= cath_lo) & (li_z < cath_hi)) | ((li_z >= anode_lo) & (li_z < anode_hi)),
            "bulk": (li_z >= bulk_zmin) & (li_z < bulk_zmax),
        }

        reference_counts = np.zeros(li_coords.shape[0], dtype=float)
        for selection in reference_sel.values():
            donor_coords = coords[selection.donor_atom_indices]
            reference_counts += donor_contact_counts(
                li_coords,
                donor_coords,
                box_x_angstrom,
                box_y_angstrom,
                use_xy_pbc,
                cutoff_sq,
            )

        for region, mask in li_region_masks.items():
            li_samples[region] += float(np.count_nonzero(mask))
            reference_cn_sum[region] += float(np.sum(reference_counts[mask]))

        for slug, selection in additive_sel.items():
            donor_coords = coords[selection.donor_atom_indices]
            additive_counts = donor_contact_counts(
                li_coords,
                donor_coords,
                box_x_angstrom,
                box_y_angstrom,
                use_xy_pbc,
                cutoff_sq,
            )
            stats = additive_stats[slug]
            for region, mask in li_region_masks.items():
                stats["li_cn_sum"][region] += float(np.sum(additive_counts[mask]))
                stats["li_contact_count"][region] += float(np.count_nonzero(additive_counts[mask] > 0.0))

    if analyzed_frames == 0:
        raise RuntimeError("No trajectory frames were analyzed")

    for stats in additive_stats.values():
        open_runs = stats["current_run_lengths"] > 0
        if np.any(open_runs):
            stats["event_lengths_frames"].extend(stats["current_run_lengths"][open_runs].tolist())
        del stats["current_run_lengths"]

    result = {
        "system_name": str(descriptor_config.get("system_name", traj_path.stem)),
        "formulation_id": str(
            args.formulation_id
            if args.formulation_id is not None
            else descriptor_config.get("formulation_id", descriptor_config.get("system_name", traj_path.stem))
        ),
        "topology_path": str(top_path),
        "trajectory_path": str(traj_path),
        "config_path": str(config_path),
        "descriptor_config_path": str(descriptor_config_path),
        "n_frames_analyzed": analyzed_frames,
        "frame_dt_ps": effective_dt_ps,
        "z_cathode_surface_angstrom": z_cathode,
        "z_anode_surface_angstrom": z_anode,
        "electrode_margin_angstrom": electrode_margin,
        "interface_width_angstrom": interface_width,
        "bulk_zmin_angstrom": bulk_zmin,
        "bulk_zmax_angstrom": bulk_zmax,
        "li_coordination_cutoff_angstrom": li_cutoff,
    }

    result.update(maybe_compute_capacitance(config, charge_log_path, args.capacitance_skip, args.capacitance_nblocks))

    reference_labels = ", ".join(selection.definition.name for selection in reference_sel.values())
    result["reference_species"] = reference_labels
    for region in ("all", "interface", "bulk"):
        result[f"reference_li_coordination_number_{region}"] = safe_ratio(reference_cn_sum[region], li_samples[region])
        result[f"li_samples_{region}"] = li_samples[region]

    interface_total_width = (cath_hi - cath_lo) + (anode_hi - anode_lo)
    for slug, stats in additive_stats.items():
        prefix = slug
        cath_density = density_from_counts(stats["cathode_count_sum"], analyzed_frames, area_nm2, cath_hi - cath_lo)
        anode_density = density_from_counts(stats["anode_count_sum"], analyzed_frames, area_nm2, anode_hi - anode_lo)
        interface_density = density_from_counts(
            stats["cathode_count_sum"] + stats["anode_count_sum"],
            analyzed_frames,
            area_nm2,
            interface_total_width,
        )
        bulk_density = density_from_counts(stats["bulk_count_sum"], analyzed_frames, area_nm2, bulk_zmax - bulk_zmin)
        residence_mean_ps, residence_p90_ps, event_count = summarize_events(
            stats["event_lengths_frames"], effective_dt_ps
        )

        result[f"{prefix}_molecule_count"] = stats["n_molecules"]
        result[f"{prefix}_cathode_number_density_1_per_nm3"] = cath_density
        result[f"{prefix}_anode_number_density_1_per_nm3"] = anode_density
        result[f"{prefix}_interface_number_density_1_per_nm3"] = interface_density
        result[f"{prefix}_bulk_number_density_1_per_nm3"] = bulk_density
        result[f"{prefix}_interface_enrichment_factor"] = safe_ratio(interface_density, bulk_density)
        result[f"{prefix}_surface_contact_probability"] = safe_ratio(
            stats["occupancy_sum"], stats["molecule_frame_count"]
        )
        result[f"{prefix}_first_layer_residence_time_ps"] = residence_mean_ps
        result[f"{prefix}_first_layer_residence_time_p90_ps"] = residence_p90_ps
        result[f"{prefix}_first_layer_event_count"] = event_count

        for region in ("all", "interface", "bulk"):
            add_cn = stats["li_cn_sum"][region]
            ref_cn = reference_cn_sum[region]
            result[f"{prefix}_li_shell_coordination_number_{region}"] = safe_ratio(add_cn, li_samples[region])
            result[f"{prefix}_li_with_contact_fraction_{region}"] = safe_ratio(
                stats["li_contact_count"][region], li_samples[region]
            )
            result[f"{prefix}_li_shell_replacement_fraction_{region}"] = safe_ratio(add_cn, add_cn + ref_cn)

        result[f"{prefix}_li_shell_replacement_interface_minus_bulk"] = (
            result[f"{prefix}_li_shell_replacement_fraction_interface"]
            - result[f"{prefix}_li_shell_replacement_fraction_bulk"]
            if np.isfinite(result[f"{prefix}_li_shell_replacement_fraction_interface"])
            and np.isfinite(result[f"{prefix}_li_shell_replacement_fraction_bulk"])
            else float("nan")
        )

    out_json_path.parent.mkdir(parents=True, exist_ok=True)
    out_csv_path.parent.mkdir(parents=True, exist_ok=True)
    out_json_path.write_text(json.dumps(result, indent=2, sort_keys=True))

    with out_csv_path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(result.keys()))
        writer.writeheader()
        writer.writerow(result)

    print(f"Frames analyzed: {analyzed_frames}")
    print(f"Descriptor row written to: {out_csv_path}")
    print(f"Descriptor JSON written to: {out_json_path}")


if __name__ == "__main__":
    main()
