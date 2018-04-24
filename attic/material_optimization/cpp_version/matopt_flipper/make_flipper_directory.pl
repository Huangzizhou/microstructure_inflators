#!/usr/bin/env perl
# Write a flipper directory to directory.js with an entry for each .js found in
# the directory.
use strict; use warnings;
use Sort::Key::Natural 'natsort';
use File::Find;
use File::Spec;

(@ARGV == 1) || die("Usage: make_flipper_directory.pl dir\n");
my $rootDir = $ARGV[0];

my $dirJS = "$rootDir/directory.js";
unlink $dirJS;

my @flippers;
find(sub { if (/\.js$/) {
               push(@flippers, File::Spec->abs2rel($File::Find::name, $rootDir));
           } }, $rootDir);

my @entries = map({ (my $name = $_) =~ tr#/#:#; $name =~ s/\.js$//; qq(\t['$name', '$_']) } natsort(@flippers));
open(DIR_OUT, '>', $dirJS);
print DIR_OUT "flippers = [\n" . join(",\n", @entries) . "\n];\n";
close(DIR_OUT);
