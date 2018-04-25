for i in stdout_*.txt; do
    python /home/fjp234/microstructures/worst_case_stress/experiments/min_autocover/extract_min_linf.py $i 0.05 0.03;
done
