#!/usr/bin/env perl
# Generate flippers (frames.js) for each .{msh,txt} pair of files in the
# specified path, writing a flipper directory to directory.js.
# The .txt file is the STDOUT of the cpp material optimization code, and the
# .msh is its output msh file.
# Usage: make_flippers_isotropic_2D.pl [path]
use strict;
use Cwd qw(abs_path cwd);
use File::Basename;
use Sort::Key::Natural 'natsort';
my $scriptDir = abs_path(dirname(__FILE__));

chdir((@ARGV >= 1) ? @ARGV[0] : cwd());

unlink 'directory.js';
my @names;
while (<*.msh>) {
    s/\.[^.]+$//;
    `$scriptDir/make_flipper.pl $_.txt $_.msh $_ > $_.js`;
    push(@names, $_);
}

my @jsFiles = map(qq(\t'$_.js'), natsort(@names));
open(DIR_OUT, '>', 'directory.js');
print DIR_OUT "flippers = [\n" . join(",\n", @jsFiles) . "\n];\n";
close(DIR_OUT);
