#!/usr/bin/env zsh
# Generate flippers (frames.js) for each .{msh,txt} pair of files in the current
# directory, writing a flipper directory to directory.js.
# The .txt file is the STDOUT of the cpp material optimization code, and the
# .msh is its output msh file.
rm -f directory.js
echo "flippers = [" > directory.js;
first=1
for i in *.msh; do
    name=${i%.msh};
    ./make_flipper.pl $name.txt $name.msh $name > $name.js;
    [[ $first == 1 ]] || echo ", " >> directory.js;
    first=0;
    echo -n "\t'$name.js'" >> directory.js
done
echo -e "\n];" >> directory.js;
