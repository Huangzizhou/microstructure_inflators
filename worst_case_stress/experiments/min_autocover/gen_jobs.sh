echo "["
expDir="$SCRATCH/wcsmin/autocover/"
for pat in {0024,0077,0646,0746,1053,1065}; do
    for num in {1..$(wc -l < tables/$pat.txt)}; do
        echo "  {'cwd': '$expDir', 'cmd': 'zsh $MICRO_DIR/worst_case_stress/experiments/min_autocover/run.sh $pat $num', 'stdout':'$expDir/stdout_${pat}_$num.txt', 'stderr': '$expDir/stderr_${pat}_$num.txt'},"
    done
done
echo "]"
