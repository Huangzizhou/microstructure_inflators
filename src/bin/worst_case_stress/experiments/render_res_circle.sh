subdiv=$1
P=$2
weight=$3

view=$MICRO_DIR/worst_case_stress/view_ptwise_measure.opt
[ $# -ge 4 ] && view=$4

suffix=""
[ $# -ge 5 ] && suffix=$5

cd /scratch/fjp234/wcs_bdry_perturb/res_circle_sweep/jvol_sweep/P$P/$subdiv
name="P${P}_jvol${weight}_sub${subdiv}"
mkdir -p $name/frames

for it in it_$name*.remeshed.msh; do
    itname=$(basename $it .remeshed.msh)
    num=${itname##*_}
    gmsh_offscreen $it $view -o $(printf "$name/frames/%04d.png" $num)
done

i=0
for img in $(lfs find $name/frames -name '*.png' | sort); do
    convert $img -gravity south \
        -stroke '#000C' -strokewidth 2 -annotate 0 "$(basename $img .png)" \
        -stroke none    -fill white    -annotate 0 "$(basename $img .png)" \
        $(printf "$name/frames/frame_%04d.png" $i)
    i=$(($i + 1))
    rm $img
done
ffmpeg -framerate 30 -i $name/frames/frame_%04d.png  -vcodec libx264 -pix_fmt yuv420p -crf 18 -y $name$suffix.mp4
rm -r $name/frames
