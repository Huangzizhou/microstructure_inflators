# usage: filter_field.sh data.txt field_no lower upper
# filter out points with field "field_no" below "lower" or above "upper"
if [[ ! $# -eq 4 ]]; then
    echo "usage: filter.sh data.txt field_no lower upper"
    exit;
fi
data=$1;
field=$2;
lower=$3;
upper=$4;

awk "((\$$field >  $lower) && (\$$field < $upper)) { print; }" $data;
