#!/usr/bin/env perl
use Parallel::ForkManager;
use List::Util qw(min max);

# Trims a collection of images consistently by cropping to the minimal bounding
# box that contains all trimmed images (with an optional margin).
(@ARGV < 1 || @ARGV > 2) && die("usage: crop.pl image_pattern [margin]\n");
my $pattern = $ARGV[0];
my $margin = (@ARGV == 2) ? $ARGV[1] : 0;

my @images = glob $pattern;
my $numImages = @images;
unlink 'bboxes.txt';
my $pm = Parallel::ForkManager->new(2);
for (my $i = 0; $i < $numImages; ++$i) {
    $pm->start and next;
    `convert ${images[$i]} -trim -format "%w %h %X %Y\n" info: >> bboxes.txt`;
    $pm->finish;
}
$pm->wait_all_children();

my (@ulx, @uly, @brx, @bry);
open(BBOX_FILE, '<', 'bboxes.txt');
while (<BBOX_FILE>) {
    my ($w, $h, $x, $y) = split(' ', $_);
    push(@ulx, $x);
    push(@uly, $y);
    push(@brx, $ulx[$i] + $w);
    push(@bry, $uly[$i] + $h);
}
close(BBOX_FILE);

# extract bb and apply margin
$bb_ulx = min(@ulx) - $margin;
$bb_uly = min(@uly) - $margin;
$bb_brx = max(@brx) + $margin;
$bb_bry = max(@bry) + $margin;

print("$bb_ulx $bb_uly $bb_brx $bb_bry\n");

my ($width, $height) = ($bb_brx - $bb_ulx, $bb_bry - $bb_uly);
for my $img (@images) {
    $pm->start and next;
    `mogrify -crop ${width}x${height}+$bb_ulx+$bb_uly +repage $img`;
    $pm->finish;
}

$pm->wait_all_children();
