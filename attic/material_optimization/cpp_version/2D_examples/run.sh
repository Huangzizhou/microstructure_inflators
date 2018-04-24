#!/usr/bin/env zsh
# Run everything!
# usage: ./run.sh [resultsDir]
resultsDir=${1-results}
if [[ -d $resultsDir ]]; then
    backupFile=${resultsDir}_old
    if [[ -e $backupFile ]]; then
        i=0;
        while [[ -e ${backupFile}_$i ]]; do
            let i++
        done
        backupFile=${backupFile}_$i
    fi
    mv $resultsDir $backupFile
fi
mkdir $resultsDir;
./generate_geometry.sh $resultsDir

# These should generate results in subdirectories "birds" and "bars"
./optimize_bars.sh $resultsDir
./optimize_birds.sh $resultsDir

# Remove input geometry
rm $resultsDir/*.obj $resultsDir/*.msh

../matopt_flipper/make_flippers.pl isotropic 2 $resultsDir
cp custom_directory.js $resultsDir/directory.js
