./LpHoleGradientComponentValidation default_job_lphole.opt 0 60 -l 0.1 -u 0.9 -a "0.0001" -m $MICRO_DIR/materials/B9Creator.material -o lphole_iterates/param_radius | tee lphole_iterates/radius.txt
./LpHoleGradientComponentValidation default_job_lphole.opt 1 60 -l 1.0 -u 5.0 -a "0.0001" -m $MICRO_DIR/materials/B9Creator.material -o lphole_iterates/param_p | tee lphole_iterates/p.txt
