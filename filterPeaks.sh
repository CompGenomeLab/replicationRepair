#!/bin/bash

awk '{if ($3-$2>1000){print}}' results/regions/R21061297-EdUrep-2hrls_combined_peaks_repdomains_org.broadPeak > results/regions/R21061297-EdUrep-2hrls_combined_peaks_repdomains_h1000.broadPeak
awk '{if ($3-$2>1000){print}}' results/regions/R21061297-EdUrep-4hrls_combined_peaks_repdomains_org.broadPeak > results/regions/R21061297-EdUrep-4hrls_combined_peaks_repdomains_h1000.broadPeak
awk '{if ($3-$2>1000){print}}' results/regions/R21061297-EdUrep-UV1-5hrls_combined_peaks_repdomains_org.broadPeak > results/regions/R21061297-EdUrep-UV1-5hrls_combined_peaks_repdomains_h1000.broadPeak
awk '{if ($3-$2>1000){print}}' results/regions/R21061297-EdUrep-UV3-5hrls_combined_peaks_repdomains_org.broadPeak > results/regions/R21061297-EdUrep-UV3-5hrls_combined_peaks_repdomains_h1000.broadPeak
awk '{if ($3-$2>1000){print}}' results/regions/R21071354-EdUrep2-2hrls2_combined_hg19_peaks_repdomains_org.broadPeak > results/regions/R21071354-EdUrep2-2hrls2_combined_hg19_peaks_repdomains_h1000.broadPeak
awk '{if ($3-$2>1000){print}}' results/regions/R21071354-EdUrep2-4hrls2_combined_hg19_peaks_repdomains_org.broadPeak > results/regions/R21071354-EdUrep2-4hrls2_combined_hg19_peaks_repdomains_h1000.broadPeak
awk '{if ($3-$2>1000){print}}' results/regions/R21071354-EdUrep2-UV1-5hrls2_combined_hg19_peaks_repdomains_org.broadPeak > results/regions/R21071354-EdUrep2-UV1-5hrls2_combined_hg19_peaks_repdomains_h1000.broadPeak
awk '{if ($3-$2>1000){print}}' results/regions/R21071354-EdUrep2-UV3-5hrls2_combined_hg19_peaks_repdomains_org.broadPeak > results/regions/R21071354-EdUrep2-UV3-5hrls2_combined_hg19_peaks_repdomains_h1000.broadPeak