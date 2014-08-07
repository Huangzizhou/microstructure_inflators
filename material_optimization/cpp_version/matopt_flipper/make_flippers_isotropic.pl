#!/usr/bin/env perl
# Generate flippers (frames.js) for each .{msh,txt} pair of files in the
# specified path, writing a flipper directory to directory.js.
# The .txt file is the STDOUT of the cpp material optimization code, and the
# .msh is its output msh file.
# Usage: make_flippers_isotropic.pl dimension [path]
use strict;
use Cwd qw(abs_path cwd);
use File::Basename;
use Sort::Key::Natural 'natsort';
my $scriptDir = abs_path(dirname(__FILE__));

(@ARGV < 1 || @ARGV > 2) && die("Usage: make_flippers_isotropic.pl dimension [directory]\n");

my $dim = $ARGV[0];
($dim ~~ [2, 3]) || die('Dimension must be 2 or 3.');

chdir((@ARGV >= 2) ? @ARGV[1] : cwd());

unlink 'directory.js';
my @names;
while (<*.msh>) {
    s/\.[^.]+$//;
    `$scriptDir/make_flipper.pl $_.txt $_.msh $_ > $_.js $dim`;
    push(@names, $_);
}

my @jsFiles = map(qq(\t'$_.js'), natsort(@names));
open(DIR_OUT, '>', 'directory.js');
print DIR_OUT "flippers = [\n" . join(",\n", @jsFiles) . "\n];\n";
close(DIR_OUT);
