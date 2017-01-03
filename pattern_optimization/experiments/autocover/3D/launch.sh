for pat in {1053,1065,0746,0646,0077,0024}; do
    pushd rounds/$pat
    python $MICRO_DIR/pattern_optimization/experiments/autocover/autocover_hpc.py 1 ../../autocover_configs/$pat.config
    popd
done
