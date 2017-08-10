cd $1
name=$(basename $1)

rm -f frame_*.png
for it in ../it_$name*.remeshed.msh; do
    itname=$(basename $it .remeshed.msh)
    num=${itname##*_}
    gmsh_offscreen $it $MICRO_DIR/worst_case_stress/view_ptwise_measure.opt -o $num.png
done

i=0
for img in $(lfs find . -maxdepth 1 -name '*.png' | sed 's/\.\///' | sort -n); do
    convert $img -gravity south \
        -stroke '#000C' -strokewidth 2 -annotate 0 "$img" \
        -stroke none    -fill white    -annotate 0 "$img" \
        $(printf "frame_%04d.png" $i)
    i=$(($i + 1))
    rm $img
done
ffmpeg -framerate 30 -i frame_%04d.png  -vcodec libx264 -pix_fmt yuv420p -crf 18 -y $name.mp4
rm *.png
