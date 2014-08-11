#!/usr/bin/env perl
# Write a flipper directory to directory.js with an entry for each .js found in
# the directory.
use strict; use warnings;
use Sort::Key::Natural 'natsort';

(@ARGV == 1) || die("Usage: make_flipper_directory.pl dir\n");
chdir($ARGV[0]);

unlink 'directory.js';
my @names = glob '*.js';

my @jsFiles = map(qq(\t'$_'), natsort(@names));
open(DIR_OUT, '>', 'directory.js');
print DIR_OUT "flippers = [\n" . join(",\n", @jsFiles) . "\n];\n";
close(DIR_OUT);
