for size in {64,128,256,512}; do
for i in input_$size.part_*; do
    part=${i#input_}
    for msh in 2D_cells/*; do
        create_pbs_noredirect.sh "$msh.$part" "cd $SCRATCH/deformed_cell_results; ./run.sh $msh $part $i" 1 4 2 0 > job_$size.$(basename $msh).$part.pbs;
    done
done
done
