#!/usr/bin/env zsh
mkdir -p init_jobs autocover_configs
ijobs=$(readlink -f init_jobs)
TDIR=$MICRO_DIR/patterns/3D/reference_wires
ACDir=$MICRO_DIR/pattern_optimization/experiments/autocover/3D
for pat in {1053,1065,0746,0646,0077,0024}; do
    patnum=$(($pat + 0))
    i=$TDIR/pattern$pat.wire
    $MICRO_DIR/pattern_optimization/GenIsosurfaceJob $i -e'1,0' -o'-0.2,0.2' > init_jobs/$pat.opt;
    cat > autocover_configs/$pat.config <<END
{
    "dim": 3,
    "pattern": $patnum,
    "material": "B9Creator",
    "jobTemplate": "$ijobs/$pat.opt",

    "targetERange": [0.15, 200],
    "targetNuRange": [-1.0, 0.5],
    "targetNSubdiv": 30,

    "numIters": 20,
    "maxSimultaneousJobs": 30,

    "mem": "16GB",
    "nprocs": 4,
    "walltime": "0:45:00",

    "args": ["--solver", "slsqp", "-V", "-O", "--TensorFitConstraint",
             "--tensor_fit_tolerance=1e-5",
             "-M", "$ACDir/coarser_3d_meshing_opts.opt",
             "--WCSWeight=0.0",
             "--JSWeight=0.0",
             "--proximityRegularizationWeight=1.0"]
}
END
done

mkdir -p rounds
for pat in {1053,1065,0746,0646,0077,0024}; do
    jobfile=init_jobs/$pat.opt
    pattern=$TDIR/pattern$pat.wire
    roundDir=rounds/$pat
    mkdir -p $roundDir
    echo "$MICRO_DIR/pattern_optimization/PatternOptimization_cli -m $MICRO_DIR/materials/B9Creator.material $jobfile -p $pattern -V --solver=levenberg_marquardt -M $ACDir/coarser_3d_meshing_opts.opt -n0 > stdout_0.txt"
    $MICRO_DIR/pattern_optimization/PatternOptimization_cli -m $MICRO_DIR/materials/B9Creator.material $jobfile -p $pattern -V --solver=levenberg_marquardt -M $ACDir/coarser_3d_meshing_opts.opt -n0 > stdout_0.txt
    if [ $? -ne 0 ]; then
        echo "Optimization failed for pattern $pat"
    fi
    python $ACDir/init_lut.py $pat $roundDir/round_0000.txt
    rm stdout_0.txt
done
