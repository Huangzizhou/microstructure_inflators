if [[ $# -ne 12 ]]; then
    echo "Usage: Usage make_jobs.sh name base_job wire material p algo solver directory procs mem hrs mins";
    exit;
fi
name=$1
baseJob=$(readlink -f $2)
wire=$(readlink -f $3)
mat=$(readlink -f $4)
p=$5
algo=$6
solver=$7
results=$(readlink -f $8)
procs=$9
mem=${10}
hrs=${11}
mins=${12}
if [[ -e $results ]]; then echo "ERROR: '$results' already exists."; exit; fi
mkdir -p $results;
if [[ ! -d $results ]]; then echo "ERROR: couldn't create results directory '$results'"; exit; fi

for (( s = 0; s <= 2; s++ )); do
    for n in {1..20}; do
        dir=$results/$name/S$s/num$n
        mkdir -p $dir
        # Lower bound init job
        create_pbs.sh "gen:$name:$algo:$s:$n:low" "$MICRO_DIR/pattern_optimization/SpacedJob $baseJob $p $n $dir/job_low_ -A $algo -S $s -p $wire -m $mat -v0" $procs $mem 1 0 > $dir/gen_low.pbs
        for (( i = 0; i < n; i++ )); do
            create_pbs.sh "run:$name:$algo:$s:$n:$i:low" "$MICRO_DIR/pattern_optimization/PatternOptimization_cli $dir/job_low_$i.opt -o $dir/job_low -A $algo -S $s -p $wire -m $mat -v0 --solver=$solver > $dir/out_low.$i.txt" $procs $mem $hrs $mins > $dir/run_low.$i.pbs
        done
        # Upper bound init job
        create_pbs.sh "gen:$name:$algo:$s:$n:up" "$MICRO_DIR/pattern_optimization/SpacedJob $baseJob $p $n $dir/job_up_ -A $algo -S $s -p $wire -m $mat -v0 -U" $procs $mem 1 0 > $dir/gen_up.pbs
        for (( i = 0; i < n; i++ )); do
            create_pbs.sh "run:$name:$algo:$s:$n:$i:up" "$MICRO_DIR/pattern_optimization/PatternOptimization_cli $dir/job_up_$i.opt -o $dir/job_up -A $algo -S $s -p $wire -m $mat -v0 --solver=$solver > $dir/out_up.$i.txt" $procs $mem $hrs $mins > $dir/run_up.$i.pbs
        done
    done
done
