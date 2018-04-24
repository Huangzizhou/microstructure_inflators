for area in {0.00064,0.00016,0.00004,0.00001}; do
    for radius in {0.25,0.5,0.75}; do
        for deg in {1,2}; do
            gnuplot -e "deg=$deg; area='$area'; radius='$radius'" lphole_wcstress_min_plot.gpi
        done
    done
done

# for f in LpHoleWCStress/*.txt; do
#     fname=${f%.txt}
#     # gnuplot -e "run='$fname';" lphole_plot_Js.gpi
#     gnuplot -e "run='$fname';" lphole_plot_wc.gpi
# done
