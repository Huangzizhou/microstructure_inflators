#!/usr/bin/env perl
# Organize runs that were all written to a single directory with colons as
# "directory separators" into an actual directory structure.
use strict; use warnings;
my $rootDir = $ARGV[0];
use File::Path qw(make_path);

chdir($rootDir);

while (<*.msh>) {
    my @components = split(':', $_);
    if (@components > 1) {
        my $dir = join('/', @components[0..$#components-1]);
        make_path($dir);
        rename($_, $dir . '/' . $components[$#components]);

        # move the .txt as well
        s/\.msh$/.txt/;
        $components[$#components] =~ s/\.msh$/.txt/;
        rename($_, $dir . '/' . $components[$#components]);
    }
}
