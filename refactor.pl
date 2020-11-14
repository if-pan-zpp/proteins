#!/usr/bin/perl

$file = 'ref/code/cg.f';

$output_file = 'cg_new.f';

open(FH1, '<', $file) or die $!;
open(FH2, '>', $output_file) or die $!;

foreach $line (<FH1>) {
  $line =~ s/kqist\((\w+),(\w+),(\w+)\)/kqist\($3,$1,$2\)/g;
  print FH2 $line;
}

close(FH1);
close(FH2);
