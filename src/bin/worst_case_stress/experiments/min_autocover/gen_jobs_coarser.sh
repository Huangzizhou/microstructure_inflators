echo "["
expDir="$SCRATCH/wcsmin/autocover_coarser_p8"
for pat in {0024,0077,0646,0746,1053,1065}; do
    for num in {1..$(wc -l < tables/$pat.txt)}; do
        echo "  {\"cwd\": \"$expDir\", \"cmd\": \"zsh /home/fjp234/microstructures/worst_case_stress/experiments/min_autocover/run_coarser_p8.sh $pat $num\", \"stdout\":\"$expDir/stdout_${pat}_$num.txt\", \"stderr\": \"$expDir/stderr_${pat}_$num.txt\"},"
    done
done
echo "]"
