echo "["
i=0
for subdiv in {64,90,128,181,256}; do
    for P in {3..12}; do
        for weight in {0,32,64,128,256}; do
            if [ $i -gt 0 ]; then 
                echo -n ", "
            fi
            echo "{\"cmd\": \"$MICRO_DIR/worst_case_stress/experiments/render_res_circle.sh $subdiv $P $weight\"}"
            i=$((i + 1))
        done
    done
done
echo "]"
