if [[ $# -ne 11 ]]; then
    echo "Usage: make_jobs.sh name base_job wire material algo solver directory procs mem hrs mins";
    exit;
fi
name=$1
baseJob=$(readlink -f $2)
wire=$(readlink -f $3)
mat=$(readlink -f $4)
algo=$5
solver=$6
results=$(readlink -f $7)
procs=$8
mem=${9}
hrs=${10}
mins=${11}
if [[ -e $results ]]; then echo "ERROR: '$results' already exists."; exit; fi
mkdir -p $results;
if [[ ! -d $results ]]; then echo "ERROR: couldn't create results directory '$results'"; exit; fi

for (( s = 0; s <= 3; s++ )); do
    for D in {0.05,0.10,0.20,0.25,0.33,0.50,1.00,16.0}; do
        dir=$results/$name/S$s/dist$D
        mkdir -p $dir;
        cmd=""
        for t in {1..20}; do
            cmd="$cmd$MICRO_DIR/pattern_optimization/RandomJob $baseJob $dir/run.$t.opt -A $algo -S $s -p $wire -m $mat -v0 -D $D"$'\n';
        done
        create_pbs.sh "gen:$name:$algo:$s:$D:$t" "$cmd" $procs $mem 1 0 >> $dir/gen.pbs
        for t in {1..20}; do
            create_pbs.sh "run:$name:$algo:$s:$D:$t" "$MICRO_DIR/pattern_optimization/PatternOptimization_cli $dir/run.$t.opt -o $dir/mesh.$t --dofOut $dir/result.$t.params -A $algo -S $s -p $wire -m $mat -v0 --solver=$solver > $dir/result.$t.txt 2>&1" $procs $mem $hrs $mins > $dir/run.$t.pbs
        done
    done
done
