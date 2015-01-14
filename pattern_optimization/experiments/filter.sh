# run after analyze.sh
# usage: filter.sh isotropy_dist
# Filters out points outside the (1 +/- isotropy) range of anisotropy measure.
if [[ ! $# -eq 2 ]]; then
    echo "usage: filter.sh isotropy_dist data.txt"
    exit;
fi
dist=$1;
data=$2;

awk "(\$4 > (1 - $dist) && \$4 < (1 + $dist)) { print; }" $data | sort -g;
