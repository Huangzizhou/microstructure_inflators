#!/usr/bin/env zsh
mkdir -p init_jobs autocover_configs
TDIR=$MICRO_DIR/patterns/2D/topologies
ACDir=$MICRO_DIR/pattern_optimization/experiments/autocover/2D
for i in $TDIR/*.obj; do
    pat=$(basename $i .obj)
    $MICRO_DIR/pattern_optimization/GenIsosurfaceJob $i -e'1,0' -o'-0.3,0.3' > init_jobs/$pat.opt;
    cat > autocover_configs/$pat.config <<END
{
    "dim": 2,
    "pattern": $(($pat + 0)),
    "material": "B9Creator",
    "jobTemplate": "$ACDir/init_jobs/$pat.opt",

    "targetERange": [0.15, 200],
    "targetNuRange": [-1.0, 1.0],
    "targetNSubdiv": 30,

    "numIters": 15,
    "maxSimultaneousJobs": 200,

    "mem": "2GB",
    "walltime": "0:15:00",

    "args": ["--solver", "levenberg_marquardt", "-V",
             "-M", "$ACDir/meshing_opts_adaptive.json",
             "--WCSWeight", "0", "--JSWeight", "1"]
}
END
done

mkdir -p rounds
for i in init_jobs/*.opt; do
    pat=$(basename $i .opt)
    roundDir=rounds/$pat
    mkdir -p $roundDir
    ../../../PatternOptimization_cli -m $MICRO_DIR/materials/B9Creator.material $i -p $TDIR/$pat.obj -V --solver=levenberg_marquardt -M meshing_opts.json -n0 > stdout_0.txt
    if [ $? -ne 0 ]; then
        echo "Optimization failed for pattern $pat"
    fi
    python init_lut.py $pat $roundDir/round_0000.txt
    rm stdout_0.txt
done
