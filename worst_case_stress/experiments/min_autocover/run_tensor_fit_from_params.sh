pat=$1
num=$2
params=$3

pattern=$MICRO_DIR/patterns/3D/reference_wires/pattern$pat.wire
db=$MICRO_DIR/worst_case_stress/experiments/min_autocover/tables/$pat.txt

line=$(head $db -n$num | tail -n1)
young=$(echo "$line" | cut -f2)
poisson=$(echo "$line" | cut -f3)
jobfile=job_${pat}_$num.opt

$PATOPT/GenIsosurfaceJob -e "$young,$poisson" -o '-0.3,0.3' -r '0.04,0.2' $pattern > $jobfile

# NOTE: no timeout on this one--we're going to run without chunking
$MICRO_DIR/worst_case_stress/WCSOptimization_cli $jobfile --JSWeight 1.0 --WCSWeight 1e-300 \
    -p $pattern -m $MICRO_DIR/materials/B9Creator.material \
    -M $MICRO_DIR/worst_case_stress/experiments/min_autocover/coarser_3d_meshing_opts.opt \
	 --TensorFitConstraint --PrintabilityConstraint \
     --params "$params" \
    -O -V -P8 -R --solver slsqp -n200
