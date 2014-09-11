#!/usr/bin/env zsh
# Run everything!
if [[ -d results ]]; then
    backupFile=results_old
    if [[ -e $backupFile ]]; then
        i=0;
        while [[ -e ${backupFile}_$i ]]; do
            let i++
        done
        backupFile=${backupFile}_$i
    fi
    mv results $backupFile
fi
mkdir results;
./generate_geometry.sh

# These should generate results in subdirectories "birds" and "bars"
./optimize_bars.sh
./optimize_birds.sh

# Remove input geometry
rm results/*.obj results/*.msh

../matopt_flipper/make_flippers.pl isotropic 2 results
cp custom_directory.js results/directory.js
