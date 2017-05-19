result_dir=$1
echo "["
i=0
for run in $(lfs find $result_dir -name '*.msh' | sed 's/it_//; s/_[0-9.][0-9.]*.remeshed.msh//' | sort -u); do
    mkdir -p $run
    if [ $i -eq 0 ]; then
        echo "{\"cmd\": \"$MICRO_DIR/worst_case_stress/experiments/postprocess_iterates.sh $run\"}"
    else
        echo ", {\"cmd\": \"$MICRO_DIR/worst_case_stress/experiments/postprocess_iterates.sh $run\"}"
    fi
    i=$((i + 1))
done
echo "]"
