#!/usr/bin/perl
#####################################
# Program: joinData.pl  -  Date: Tue Apr 14 15:15:49 BRT 2015
# Autor: Rubens Magalhaes
#
# License: GPL - http://www.gnu.org/licenses/gpl.html
#
#####################################
open (IN0,"<$ARGV[0]");
open (IN1,"<$ARGV[1]");
while (<IN0>) {
    chomp $_;
    @tmp = split (/\t+/, $_);
    @tmp1 = split (/\;/, $tmp[8]);
    $genes{"$tmp[3] $tmp[4] $tmp[6]"} = "$tmp1[1] $tmp1[2]";
}

while (<IN1>) {
    chomp $_;
    @tmp = split (/\s+/, $_);
    print "$_ ##";
    $tmp[10] =~ s/"//g;
    print $genes{"$tmp[8] $tmp[9] $tmp[10]"};
    print "\n";
}
