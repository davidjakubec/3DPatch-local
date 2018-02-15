#!/usr/bin/perl

use strict;
use warnings;
use Bio::HMM::Logo;

my $hmmfile = shift;

my $logo = Bio::HMM::Logo->new({hmmfile => $hmmfile});
my $json = $logo->as_json;

open(F, ">", "logo.json");
print F $json;
close(F);
