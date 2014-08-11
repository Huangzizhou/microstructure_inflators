#!/usr/bin/env perl
# Generate flippers (frames.js) for each .{msh,txt} pair of files in the
# specified path.
# The .txt file is the STDOUT of the cpp material optimization code, and the
# .msh is its output msh file.
# Usage: make_flippers.pl material dimension [path [filter]]
use strict; use warnings;
use Cwd qw(abs_path cwd);
use File::Basename;
my $scriptDir = abs_path(dirname(__FILE__));

(@ARGV < 2 || @ARGV > 4) && die("Usage: make_flippers.pl iso|orthotropic dimension [directory [filter]]\n");

my $material = $ARGV[0];
my $dim = $ARGV[1];

($material ~~ ['isotropic', 'orthotropic']) || die('Material must be isotropic or orthotropic');
($dim ~~ [2, 3]) || die('Dimension must be 2 or 3.');

chdir((@ARGV >= 3) ? $ARGV[2] : cwd());

my $properties;
if ($material eq 'isotropic') {
       $properties = 'E nu';
}
else { $properties = ($dim == 2) ? 'E_x E_y nu_yx mu'
                                 : 'E_x E_y E_z nu_yx nu_zx nu_zy mu_yz mu_zx mu_xy';
}
my $commaProperties = $properties;
$commaProperties =~ s/ /, /g;

unlink 'directory.js';
my @names = glob (((@ARGV >= 4) ? $ARGV[3] : '*') . '.msh');
@names = map {local $_ = $_; s/\.[^.]+$//; $_} @names;
for (my $i = 0; $i < @names; ++$i) {
    my $name = $names[$i];
    `$scriptDir/make_flipper.pl $name.txt $name.msh $name $dim '$properties' '$commaProperties' > $name.js`;
    print "$name\t(" . ($i + 1) . '/' . @names . ")\n";
}
