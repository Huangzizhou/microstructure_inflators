#!/usr/bin/env perl
# Generate flippers (frames.js) for each .{msh,txt} pair of files in the
# specified path.
# The .txt file is the STDOUT of the cpp material optimization code, and the
# .msh is its output msh file.
# Usage: make_flippers.pl material dimension [path [filter]]
use strict; use warnings;
use Cwd qw(abs_path cwd);
use File::Basename;
use File::Find;
use File::chdir;
use File::Spec;
my $scriptDir = abs_path(dirname(__FILE__));

(@ARGV < 2 || @ARGV > 4) && die("Usage: make_flippers.pl iso|orthotropic dimension [directory [filter]]\n");

my $material = $ARGV[0];
my $dim = $ARGV[1];

($material ~~ ['isotropic', 'orthotropic']) || die('Material must be isotropic or orthotropic');
($dim ~~ [2, 3]) || die('Dimension must be 2 or 3.');

my $rootDir = (@ARGV >= 3) ? $ARGV[2] : cwd();

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
my @mshFiles;
my $filter = ((@ARGV >= 4) ? $ARGV[3] : '') . '.*\.msh';
find(sub { my $path = File::Spec->abs2rel($File::Find::name, $rootDir);
           if ($path =~ /$filter/) { push(@mshFiles, $path); } }, $rootDir);

for (my $i = 0; $i < @mshFiles; ++$i) {
    my ($name, $dir, $suffix) = fileparse($mshFiles[$i], qr/\.[^.]*/);
    local $CWD = "$rootDir/$dir";
    (my $runName = $dir) =~ tr#/#:#;
    $runName .= $name;
    `$scriptDir/make_flipper.pl $name.txt $name.msh '$runName' $dim '$properties' '$commaProperties' > $name.js`;
    print "$runName\t(" . ($i + 1) . '/' . @mshFiles . ")\n";
}
