pat=$1
num=$2

pattern=$MICRO_DIR/patterns/3D/reference_wires/pattern$pat.wire
db=$MICRO_DIR/worst_case_stress/experiments/min_autocover/tables/$pat.txt

line=$(head $db -n$num | tail -n1)
young=$(echo "$line" | cut -f2)
poisson=$(echo "$line" | cut -f3)
params=$(echo "$line" | cut -f6-)
jobfile=job_${pat}_$num.opt

$PATOPT/GenIsosurfaceJob -e "$young,$poisson" -o '-0.3,0.3' -r '0.04,0.2' $pattern > $jobfile

# We are running many of these jobs with chunking, so we must kill each job
# before it takes 1/Nth of the requested time otherwise the later jobs might be killed.
echo "$MICRO_DIR/worst_case_stress/WCSOptimization_cli $jobfile --JSWeight 0.0 --WCSWeight 1.0 \
    -p $pattern -m $MICRO_DIR/materials/B9Creator.material \
    -M $MICRO_DIR/worst_case_stress/experiments/min_autocover/finer_3d_meshing_opts.opt \
	 --TensorFitConstraint --PrintabilityConstraint \
     --params "$params" \
    -O -V -P12 -R --solver slsqp -n200"
