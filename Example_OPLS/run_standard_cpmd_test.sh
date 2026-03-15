#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
WORKFLOW_DIR="$SCRIPT_DIR/LiPF6_EC_DMC_thick_Li_electrode/three_stage_cp_workflow"

MODE="${1:-smoke}"
PLATFORM="${2:-CPU}"

cd "$WORKFLOW_DIR"

if [[ "$MODE" == "smoke" ]]; then
  python run_three_stage_cp_workflow.py \
    --platform "$PLATFORM" \
    --neutral-relax-steps 5000 \
    --cp-activation-steps 5000 \
    --production-steps 10000 \
    --report-interval 500
elif [[ "$MODE" == "standard" ]]; then
  python run_three_stage_cp_workflow.py \
    --platform "$PLATFORM" \
    --neutral-relax-steps 200000 \
    --cp-activation-steps 100000 \
    --production-steps 500000 \
    --report-interval 1000
else
  echo "Usage: bash run_standard_cpmd_test.sh [smoke|standard] [CPU|CUDA|OpenCL|Reference]" >&2
  exit 1
fi

python analyze_electrode_total_charge.py \
  --input outputs/electrode_total_charge_timeseries.dat \
  --output outputs/electrode_total_charge_summary.txt

echo "Finished $MODE run in $WORKFLOW_DIR"
