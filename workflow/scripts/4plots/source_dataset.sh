#!/bin/bash

####Choose the data####

#data_name=( chromhmm.repdomains.hela_windows.bed chromhmm.repdomains.hela_windows.bed genome.hg19_windows_100kb.bed genome.hg19_windows_1Mb.bed genome.hg19_windows_20kb.bed hi.rfd.hela_windows_201_100.bed hi.rfd.hela_windows_201_1000.bed hi.rfd.hela_windows_201_10000.bed initiation.zones.hela_windows_201_100.bed initiation.zones.hela_windows_201_1000.bed initiation.zones.hela_windows_201_10000.bed repdomains.hela.bed repdomains.hela_windows_201_100.bed repdomains.hela_windows_201_1000.bed repdomains.hela_windows_201_10000.bed s1.to.s2.20kb.hela_windows_201_100.bed s1.to.s2.20kb.hela_windows_201_1000.bed sns.seq.hela.rep.1.2_windows_201_100.bed sns.seq.hela.rep.1.2_windows_201_1000.bed initiation.zones.repdomains.hela_windows_201_100.bed hi.rfd.repdomains.hela_windows_201_100.bed s1.to.s2.20kb.repdomains.hela_windows_201_100.bed sns.seq.repdomains.hela.rep.1.2_windows_201_100.bed )

#dataset=( chromhmm_windows chromhmm_windows_chr genome_100kb genome_1Mb genome_20kb hiRFD_windows_201_100 hiRFD_windows_201_1000 hiRFD_windows_201_10000 inZones_windows_201_100 inZones_windows_201_1000 inZones_windows_201_10000 repdomains repdomains_windows_201_100 repdomains_windows_201_1000 repdomains_windows_201_10000 S1toS2_windows_201_100 S1toS2_windows_201_1000 sns_seq_windows_201_100 sns_seq_windows_201_1000 inZones_repdomains_windows_201_100 hiRFD_repdomains_windows_201_100 S1toS2_repdomains_windows_201_100 sns_seq_repdomains_windows_201_100 )

#combine_options=( "-strand T" "-strand T -chr T" "" "" "" "-strand T" "-strand T" "-strand T" "" "" "" "" "" "" "" "-strand T" "-strand T" "" "" "" "-strand T" "-strand T" "" ) # -score T -strand T -chr T 

#### alternatives ####

data_name=( )

dataset=( )

combine_options=( )

##### 
