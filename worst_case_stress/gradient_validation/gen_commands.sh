for deg in {1,2}; do
    for pnorm in {10..60}; do
        for area in {0.00001,0.00004,0.00016,0.00064}; do
            for radius in {0.25,0.5,0.75}; do
                echo "\$MICRO_DIR/worst_case_stress/LpHoleGradientComponentValidation <(\$MICRO_DIR/worst_case_stress/lphole_radius_job.sh $radius) 1 240 -l 0.25 -u 2.5 -d$deg -P $(($pnorm / 10.0)) -a $area -m \$MICRO_DIR/materials/B9Creator.material >! d${deg}_p${pnorm}_a${area}_r$radius.txt"
            done
        done
    done
done
