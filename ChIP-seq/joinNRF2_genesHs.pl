#!/usr/bin/perl
#####################################
# Program: joinNRF2_genesHs.pl  -  Date: Mon Sep  1 15:49:10 BRT 2014
# Autor: Rubens Magalhaes
#
# License: GPL - http://www.gnu.org/licenses/gpl.html
#
#####################################
open (IN0,"<$ARGV[0]");
open (IN1,"<$ARGV[1]");
while (<IN0>) {
    chomp $_;
    @tmp = split (/\s+/, $_);
    if ($tmp[6] eq "+") {
	$chromP{$tmp[0]} .= "$tmp[3]###$_|";
    }
    elsif ($tmp[6] eq "-") {
	$chromM{$tmp[0]} .= "$tmp[4]###$_|";
    }
}
while (<IN1>) {
    chomp $_;
    @tmp = split (/\s+/, $_);
    $chr = $tmp[0];
    $start = $tmp[1];
    $end = $tmp[2];
    $resP = get_distPlus($start,$end, $chr, $_);
    $resM = get_distMinus($start,$end, $chr, $_);
    if ($resP == 1 || $resM == 1) {
	print "-"x200;
	print "\n";
    }
}
sub get_distPlus {
    local $tf_start = $_[0] - 600;
    local $tf_end = $_[1] + 600;
    local $tf_chr = $_[2];
    local $tf_all = $_[3];
    local $print = 0;
    local $avg_pos = $tf_start + ($tf_end - $tf_start)/2;
    if($chromP{$tf_chr} ne "") {
	@coordChr = split (/\|/, $chromP{$tmp[0]});
	@coordChrP = sort {$a <=> $b} @coordChr;
	#@coordChrP = sort @coordChr;
	while (@coordChrP) {
	    $tmp = shift (@coordChrP);
	    @posit = split (/\#\#\#/, $tmp);
	    for ($i = $tf_start; $i <= $tf_end; $i++) {
		if ($posit[0] == $i) {
		    $dist = $i - $avg_pos;
		    print "$tf_all |+ dist=> $dist| $posit[1]\n";		    
		    $print = 1;
		}
	    }
	}
    }
    return $print;
}
sub get_distMinus {
    local $tf_start = $_[0] - 600;
    local $tf_end = $_[1] + 600;
    local $tf_chr = $_[2];
    local $tf_all = $_[3];
    local $print = 0;
    local $avg_pos = $tf_start + ($tf_end - $tf_start)/2;
    if($chromM{$tf_chr} ne "") {
	@coordChr = split (/\|/, $chromM{$tmp[0]});
	@coordChrM = sort {$a <=> $b} @coordChr;
	while (@coordChrM) {
	    $tmp = shift (@coordChrM);
	    @posit = split (/\#\#\#/, $tmp);
	    for ($i = $tf_start; $i <= $tf_end; $i++) {
		if ($posit[0] == $i) {
		    $dist = $avg_pos - $i;
		    print "$tf_all |- dist=> $dist| $posit[1]\n";		    
		    $print = 1;
		}
	    }
	}
    }
    return $print;
}
