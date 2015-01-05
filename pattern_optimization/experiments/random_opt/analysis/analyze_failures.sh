# To be run from a run's root results directory
# (e.g. $SCRATCH/pattern_optimization/random/brick5.01_02_15.05:23)
mkdir -p analysis
find . -name run\*.pbs > analysis/all_jobs.txt
truncate --size 0 analysis/gen_{success,failure}.txt
for i in $(find . -name run\*.pbs); do
    if  [[ -e ${i%pbs}opt ]]; then
        echo $i >> analysis/gen_success.txt
    else
        echo $i >> analysis/gen_failure.txt
    fi
done
tooShortGenErrorFile=$(mktemp)
tetgenBugGenErrorFile=$(mktemp)
grep "to ensure inflation" -A1 gen* | grep -v "^--$" | paste - - > $tooShortGenErrorFile
grep "Tetgen Error" -A5 gen* | grep -v "^--$" | paste - - - - > $tetgenBugGenErrorFile
truncate --size 0 analysis/gen_failure_{too_short,tetgen_bug,unrecognized}.txt
for i in $(cat analysis/gen_failure.txt); do
    if test -n "$(grep ${i%pbs}opt $tooShortGenErrorFile)"; then
        echo $i >> analysis/gen_failure_too_short.txt
    else
        if test -n "$(grep ${i%pbs}opt $tetgenBugGenErrorFile)"; then
            echo $i >> analysis/gen_failure_tetgen_bug.txt
        else
            echo "UNRECOGNIZED FAILURE: $i"
            echo $i >> analysis/gen_failure_unrecognized.txt
        fi
    fi
done

truncate --size 0 analysis/run_{success,failure}.txt
truncate --size 0 analysis/run_failure_{too_short,tetgen_bug,min_bdry_loop,inflator_segfault,empty,exceeded_time,unrecognized}.txt
for i in $(cat analysis/gen_success.txt); do
    resultFile=$(echo $i | sed 's/run.\([0-9]*\).pbs/result.\1.txt/')
    if test -n "$(tail $resultFile -n1 | grep 'Ceres')"; then
        echo $i >> analysis/run_success.txt;
        continue;
    fi;
    echo $i >> analysis/run_failure.txt

    if grep -q "Edge .* is too short" $resultFile; then
        echo $i >> analysis/run_failure_too_short.txt;
        continue;
    fi;

    if grep -q "Tetgen internal bug" $resultFile; then
        echo $i >> analysis/run_failure_tetgen_bug.txt
        continue;
    fi;

    if grep -q "Unable to match min boundary loops" $resultFile; then
        echo $i >> analysis/run_failure_min_bdry_loop.txt
        continue;
    fi;

    if grep -q "Empty inflated geometry" $resultFile; then
        echo $i >> analysis/run_failure_empty.txt
        continue;
    fi;

    # Note: this only works if the annoying extra status reporting prints are enabled.
    if grep -q "^Inflating$" <(tail -n1 $resultFile); then
        if grep $resultFile *.e[0-9]* | grep -q "Segmentation fault"; then
            echo $i >> analysis/run_failure_inflator_segfault.txt
            continue;
        fi
    fi;

    # Assume quits during homogenizing are due to exceeding time (so far this has been the case)
    # Note: this only works if the annoying extra status reporting prints are enabled.
    if grep -q "^Homogenizing$" <(tail -n1 $resultFile); then
        echo $i >> analysis/run_failure_exceeded_time.txt
        continue;
    fi;

    echo "UNRECOGNIZED FAILURE: $i"
    echo $i >> analysis/run_failure_unrecognized.txt
done

numJobs=$(wc -l < analysis/all_jobs.txt);
genSuccesses=$(wc -l < analysis/gen_success.txt);
echo -e "Job Generation:\t$genSuccesses/$numJobs\t($(bc <<< "scale=2; 100 * ($numJobs - $genSuccesses)/$numJobs")% failure)"
for i in {too_short,tetgen_bug,unrecognized}; do
    count=$(wc -l < analysis/gen_failure_$i.txt)
    echo -e "\terror type $i:\t$count\t($(bc <<< "scale=2; 100 * $count / $numJobs")%)"
done

runSuccesses=$(wc -l < analysis/run_success.txt);
echo -e "Optimization:\t$runSuccesses/$genSuccesses\t($(bc <<< "scale=2; 100 * ($genSuccesses - $runSuccesses)/$genSuccesses")% failure)"
for i in {too_short,tetgen_bug,min_bdry_loop,inflator_segfault,empty,exceeded_time,unrecognized}; do
    count=$(wc -l < analysis/run_failure_$i.txt)
    echo -e "\terror type $i:\t$count\t($(bc <<< "scale=2; 100 * $count / $genSuccesses")%)"
done
