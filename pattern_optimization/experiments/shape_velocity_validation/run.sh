mkdir -p results
for p in {0..8}; do ./ShapeVelocityValidation -p ../../../patterns/3D/reference_wires/pattern0746.wire -l0.03 -u0.1 ../../default_job_isoinflator.opt $p -n 20 -IV 2>&1 -o  direct_nsv_$p | tee results/direct_nsv_$p.txt; done
for p in {0..8}; do ./ShapeVelocityValidation -p ../../../patterns/3D/reference_wires/pattern0746.wire -l0.03 -u0.1 ../../default_job_isoinflator.opt $p    20 -IV 2>&1 -o project_nsv_$p | tee results/projected_nsv_$p.txt; done
cd results
for prefix in {direct_nsv_,projected_nsv_}; do for i in {0..8}; do gnuplot -e "prefix='$prefix'; paramNum=$i" plot.gpi; done; done
