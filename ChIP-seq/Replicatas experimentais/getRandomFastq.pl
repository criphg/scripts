#!/usr/bin/perl
#####################################
# Program: getRandomFastq.pl  -  Date: Qua Mar  2 13:31:00 BRT 2016
# Autor: Pedro A. F. Galante
#
# License: GPL - http://www.gnu.org/licenses/gpl.html
#
#####################################

#perl programa <infileR1> <infileR2> <outfileR1> <outfileR2>

open (IN0,"<$ARGV[0]"); # R1
open (IN1,"<$ARGV[1]"); # R2
$countR1 = 0;
$countR2 = 0;
while (<IN0>) {
    chomp $_;
    $seq = "";
    if (/^\@NS/) {
	$header = $_;
	$seq .= <IN0>;
	$seq .= <IN0>;
	$seq .= <IN0>;
	$dataR1[$countR1] = "$header\n$seq";
	$countR1++;
    }
}
while (<IN1>) {
    chomp $_;
    $seq = "";
    if (/^\@NS/) {
	$header = $_;
	$seq .= <IN1>;
	$seq .= <IN1>;
	$seq .= <IN1>;
	$dataR2[$countR2] = "$header\n$seq";
	$countR2++;
    }
}

#modificado para atribuir o nome dos arquivos output
open (OUT1, ">$ARGV[2]");
open (OUT2, ">$ARGV[3]");
for ($i = 0; $i < 1070000; $i++) {
    $rand = int(rand($countR1));
    if ($selected{$rand} eq "") {	
	print OUT1 $dataR1[$rand];
	print OUT2 $dataR2[$rand];
	$selected{$rand} = 1;
    }
    else {
	$i--;
    }
}
