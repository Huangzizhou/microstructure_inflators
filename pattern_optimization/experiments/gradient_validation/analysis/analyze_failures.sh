# To be run from a run's root results directory
# (e.g. $SCRATCH/pattern_optimization/random/brick5.01_02_15.05:23)
mkdir -p analysis
truncate --size 0 analysis/all_gen.txt
find . -name gen.pbs > analysis/all_gen.txt
truncate --size 0 analysis/gen_{success,failure}.txt
truncate --size 0 analysis/gen_failure_unrecognized.txt
for i in $(find . -name gen.pbs); do
    if [[ -e ${i%gen.pbs}run_0.opt ]]; then
        echo $i >> analysis/gen_success.txt
    else
        echo $i >> analysis/gen_failure.txt
    fi;
done

for i in $(cat analysis/gen_failure.txt); do
    echo "UNRECOGNIZED FAILURE: $i"
done

truncate --size 0 analysis/run_{success,failure_{success,tetgen_bug,exceeded_time,inflator_segfault,unrecognized}}.txt
totalSweepJobs=0
for i in $(cat analysis/gen_success.txt); do
    for r in ${i%gen.pbs}run.*.pbs; do
        totalSweepJobs=$(( $totalSweepJobs + 1 ));
        res=$(sed 's/run/result/; s/pbs/txt/' <<< $r);
        if tail -n1 $res | grep -q "Full time"; then
            echo $r >> analysis/run_success.txt
            continue
        fi

        if grep -q "Tetgen internal bug" $res; then
            echo $r >> analysis/run_failure_tetgen_bug.txt
            continue
        fi

        # Assume quits during homogenizing are due to exceeding time (so far this has been the case)
        # Note: this only works if the annoying extra status reporting prints are enabled.
        if tail -n1 $res | grep -q "Homogenizing"; then
            echo $r >> analysis/run_failure_exceeded_time.txt
            continue
        fi

        # Note: this only works if the annoying extra status reporting prints are enabled.
        if grep -q "^Inflating$" <(tail -n1 $res); then
            if grep $res *.e[0-9]* | grep -q "Segmentation fault"; then
                echo $r >> analysis/run_failure_inflator_segfault.txt
                continue;
            fi
        fi;

        echo "UNRECOGNIZED FAILURE: $r"
        echo $r >> analysis/run_failure_unrecognized.txt
    done
done

numGenJobs=$(wc -l < analysis/all_gen.txt);
genSuccesses=$(wc -l < analysis/gen_success.txt);
echo -e "Job Generation:\t$genSuccesses/$numGenJobs\t($(bc <<< "scale=2; 100 * ($numGenJobs - $genSuccesses)/$numGenJobs")% failure)"
for i in unrecognized; do
    count=$(wc -l < analysis/gen_failure_$i.txt)
    echo -e "\terror type $i:\t$count\t($(bc <<< "scale=2; 100 * $count / $numGenJobs")%)"
done

runSuccesses=$(wc -l < analysis/run_success.txt);
echo -e "Sweep:\t$runSuccesses/$totalSweepJobs\t($(bc <<< "scale=2; 100 * ($totalSweepJobs - $runSuccesses)/$totalSweepJobs")% failure)"
for i in {tetgen_bug,exceeded_time,inflator_segfault,unrecognized}; do
    count=$(wc -l < analysis/run_failure_$i.txt)
    echo -e "\terror type $i:\t$count\t($(bc <<< "scale=2; 100 * $count / $totalSweepJobs")%)"
done
