#!/usr/bin/env zsh
mkdir -p init_jobs autocover_configs
ijobs=$(readlink -f init_jobs)
TDIR=$MICRO_DIR/patterns/2D/topologies
ACDir=$MICRO_DIR/pattern_optimization/experiments/autocover/2D
for i in $TDIR/*.obj; do
    pat=$(basename $i .obj)
    $MICRO_DIR/pattern_optimization/GenIsosurfaceJob $i -e'1,0' -r'0.02,0.2' -o'-0.3,0.3' > init_jobs/$pat.opt;
    # $MICRO_DIR/pattern_optimization/GenIsosurfaceJob $i -e'1,0' -r'0.04,0.2' -o'-0.3,0.3' > init_jobs/$pat.opt;
    cat > autocover_configs/$pat.config <<END
{
    "dim": 2,
    "pattern": $(($pat + 0)),
    "material": "B9Creator",
    "jobTemplate": "$ijobs/$pat.opt",

    "targetERange": [0.15, 200],
    "targetNuRange": [-1.0, 1.0],
    "targetNSubdiv": 30,

    "numIters": 20,
    "maxSimultaneousJobs": 50,
    "singleClosestInit": true,

    "mem": "2GB",
    "walltime": "0:10:00",

    "args": ["--solver", "slsqp", "-V", "-O", "--TensorFitConstraint",
             "--tensor_fit_tolerance=1e-5",
             "-M", "$ACDir/meshing_opts_adaptive.json",
             "--WCSWeight=0.0",
             "--JSWeight=0.0",
             "--proximityRegularizationWeight=1.0"]
}
END
done

mkdir -p rounds
for i in init_jobs/*.opt; do
    pat=$(basename $i .opt)
    roundDir=rounds/$pat
    mkdir -p $roundDir
    echo "$MICRO_DIR/pattern_optimization/PatternOptimization_cli -m $MICRO_DIR/materials/B9Creator.material $i -p $TDIR/$pat.obj -V --solver=levenberg_marquardt -M $ACDir/meshing_opts.json -n0 > stdout_0.txt"
    $MICRO_DIR/pattern_optimization/PatternOptimization_cli -m $MICRO_DIR/materials/B9Creator.material $i -p $TDIR/$pat.obj -V --solver=levenberg_marquardt -M $ACDir/meshing_opts.json -n0 > stdout_0.txt
    if [ $? -ne 0 ]; then
        echo "Optimization failed for pattern $pat"
    fi
    python $ACDir/init_lut.py $pat $roundDir/round_0000.txt
    rm stdout_0.txt
done
