import numpy as np
pa = [0.6672992017706555,      0.7332146805779556,      0.4110261537397022, 0.3567761454403041,      0.02643254526871755,     0.1999571975522177, 0.1999649950360378,      0.01501359305328773,     0.0266869303815188, 0.029567810355904,       0.01,    0.01,    0.08299817740886838, 0.08720194210402653]
pa = np.array(pa)
pb = [0.7999999999998425,      0.7999999999999974,      0.3125550670520073, 0.3315163624067535,      0.02353853774294269,     0.1999999999999872, 0.1999999999999992,      0.01500000000000203,     0.01569595803363666, 0.05725613693259792,     0.01,    0.007812500000029914, 0.09977315880265886, 0.097265831915867]
pb = np.array(pb)
lerp = lambda f : pa + (pb - pa) * f
params = [lerp(f) for f in np.arange(0.0, 1.0, 0.05)]

for i, p in enumerate(params):
    print '$MICRO_DIR/worst_case_stress/WCSOptimization_cli --JSWeight 0 --WCSWeight 0 -m $MICRO_DIR/materials/B9Creator.material -p $MICRO_DIR/Luigi/wireinflator2D/meshes/octa_cell.obj <(./target_isotropic.sh 0.288 0.0494)  -M <($MICRO_DIR/worst_case_stress/gradient_validation/octacell_iso/mesh_options.sh 1e-4 2048) -V -d2 --solver=slsqp -n100 -CO --proximityRegularizationWeight 1.0 --proximityRegularizationTarget "{}" --params "{}" | tee step_{}.txt'.format(
            " ".join(map(str, p)), " ".join(map(str, params[0])), i)
