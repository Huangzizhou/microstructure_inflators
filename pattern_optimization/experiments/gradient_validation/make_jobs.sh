if [[ $# -ne 13 ]]; then
    echo "Usage: make_jobs.sh name base_job wire material algo comp_idx lower_bd upper_bd directory procs mem hrs mins";
    exit;
fi
name=$1
baseJob=$(readlink -f $2)
wire=$(readlink -f $3)
mat=$(readlink -f $4)
algo=$5
comp_idx=$6
lower_bd=$7
upper_bd=$8
results=$(readlink -f $9)
procs=${10}
mem=${11}
hrs=${12}
mins=${13}

# Allow writing to existing directory so different ranges can be run in the same directory.
# if [[ -e $results ]]; then echo "ERROR: '$results' already exists."; exit; fi
mkdir -p $results;
if [[ ! -d $results ]]; then echo "ERROR: couldn't create results directory '$results'"; exit; fi

range="comp$comp_idx,$lower_bd,$upper_bd"
for (( s = 0; s <= 3; s++ )); do
    for vol in {'0.001','0.0005','0.00025','0.000125','0.0000625'}; do
        dir="$results/$name/$algo/S$s/vol$vol";
        mkdir -p $dir;
        # Generate the target at the middle of the full parameter range (this is independent of the sweep range passsed to this script)
        cmd="$MICRO_DIR/pattern_optimization/SpacedJob $baseJob $comp_idx 1 $dir/run_ -A $algo -S $s -p $wire -m $mat -v$vol";
        create_pbs.sh "gen:$name:$algo:$s:$vol" "$cmd" $procs $mem 0 30 > $dir/gen.pbs;
        # positional args: job.opt component_idx -l lower_bd -u upper_bd nsamples
        cmd="$MICRO_DIR/pattern_optimization/GradientComponentValidation $dir/run_0.opt $comp_idx -l$lower_bd -u$upper_bd 801 --dofOut $dir/$range.params -A $algo -S $s -p $wire -m $mat -v$vol > $dir/result.$range.txt 2>&1";
        create_pbs.sh "run:$name:$algo:$s:$vol:$range" "$cmd" $procs $mem $hrs $mins > $dir/run.$range.pbs;
    done
done
