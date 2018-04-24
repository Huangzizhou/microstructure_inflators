echo -ne "sample\t"
echo -ne "init_young\tinit_poisson\topt_young\topt_poisson\t"
echo -ne "init_wcs\topt_wcs\t"
echo -ne "init_zcompress_min\tinit_zcompress_max\topt_zcompress_min\topt_zcompress_max\t"
echo -ne "init_zcompress_sim_min\tinit_zcompress_sim_max\topt_zcompress_sim_min\topt_zcompress_sim_max\t"
echo -e  "init_zcompress_sim_2x2x1_min\tinit_zcompress_sim_2x2x1_max\topt_zcompress_sim_2x2x1_min\topt_zcompress_sim_2x2x1_max"

directories=([0-9]*)
if [ $# -ne 0 ]; then
    directories=("$@")
fi

for dir in $directories; do
    echo -ne "$dir\t"
  ( grep 'Young\|v_yx' $dir/init.homog.txt | cut -f2
    grep 'Young\|v_yx' $dir/opt.homog.txt  | cut -f2
    $MeshFEM/tools/msh_processor $dir/init.wcs.msh -e "Pointwise WCS" --max
    $MeshFEM/tools/msh_processor  $dir/opt.wcs.msh -e "Pointwise WCS" --max
    $MeshFEM/tools/msh_processor $dir/init.zcompress.msh           -e "stress" -l --dup --max --max --reverse --min --min --applyAll --print
    $MeshFEM/tools/msh_processor  $dir/opt.zcompress.msh           -e "stress" -l --dup --max --max --reverse --min --min --applyAll --print
    $MeshFEM/tools/msh_processor $dir/init.zcompress.sim.msh       -e "stress" -l --dup --max --max --reverse --min --min --applyAll --print
    $MeshFEM/tools/msh_processor  $dir/opt.zcompress.sim.msh       -e "stress" -l --dup --max --max --reverse --min --min --applyAll --print
    $MeshFEM/tools/msh_processor $dir/init.zcompress.sim_2x2x1.msh -e "stress" -l --dup --max --max --reverse --min --min --applyAll --print
    $MeshFEM/tools/msh_processor  $dir/opt.zcompress.sim_2x2x1.msh -e "stress" -l --dup --max --max --reverse --min --min --applyAll --print
    ) | tr '\n' '\t'
    echo ""
done
