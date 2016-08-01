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

	$tmp1[2] =~ s/\+/ /g;
	$tmp1[2] =~ s/%2C//g;
	$tmp1[2] =~ s/%29//g;
	$tmp1[2] =~ s/%28//g;
	$tmp1[1] =~ s/Name=//g;
	$tmp1[2] =~ s/description=//g;	
    	$genes{"$tmp[3] $tmp[4] $tmp[6]"} = "\t$tmp1[1]\t$tmp1[2]";
}
$a = "feature_strand";
$b = "insideFeature";
$c = "distancetoFeature";
$genes{"$a $b $c"} = "\tID\tName"; 


while (<IN1>) {
    chomp $_;
    @tmp = split (/\s+/, $_);
    print "$_";
    $tmp[12] =~ s/"//g;
	$tmp[10] =~ s/"//g;
	$tmp[11] =~ s/"//g;

    print $genes{"$tmp[10] $tmp[11] $tmp[12]"};
    print "\n";

}
