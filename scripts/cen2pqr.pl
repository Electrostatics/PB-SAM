#!/usr/bin/perl -w
use strict;
my $decimal = '-?\d+\.?\d*';

my $atom = 1; 
my $atomname = "O";
my $aaname = "CEN";
my $aa;
my $x, my $y, my $z, my $r;

open(READ, "$ARGV[0]") || die;

while(my $line = <READ>)
{

    ($x,$y,$z, $r) = split(/\s+/, $line);
   
    $aa = $atom;

    printf "%-5s  %4d  %-4s %3s   %3d    %8.3f %8.3f %8.3f %5.2f %8.4f\n",
           "ATOM",$atom, $atomname,$aaname, $aa, $x,  $y,  $z,  0.0 ,  $r;

    $atom++;

}
