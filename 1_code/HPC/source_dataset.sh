#!/bin/bash

####Data Info####

#dataset=(repdomain EDTZ_EUTZ ERD_LRD_windows HiRFD_windows InZones_windows rep1_SNSseq rep2_SNSseq rep1_SNSseq_windows rep2_SNSseq_windows S1toS2 S1toS2_windows ChrStates_comb_windows ChrStates_segway_windows ChrStates_chromhmm_windows)

#data_name=(GSE53984_GSM923449_Helas3_Rep1_segments_updated.bed GSE53984_GSM923449_Helas3_Rep1_segments_EDTZ_EUTZ_sort.bed GSE53984_GSM923449_Helas3_Rep1_segments_ERD_LRD_windows.bed 20160224_Segment_HMM_Flat_15kb_Hela_Combined_V6_Raw_LT0_1_Selected_HighRFD_windowed_full.bed InitiationZones_Hela_updated_windows.bed Besnard_HeLa_rep1_peaks_updated.bed Besnard_HeLa_rep2_peaks_updated.bed Besnard_HeLa_rep1_just_peaks_updated_windows.bed Besnard_HeLa_rep2_just_peaks_updated_windows.bed S1toS2_peak_100kb.bed S1toS2_exact_match_100kb_windowed.bed wgEncodeAwgSegmentationCombinedHelas3_repdomains.bed wgEncodeAwgSegmentationSegwayHelas3_repdomains.bed wgEncodeAwgSegmentationChromhmmHelas3_repdomains.bed)

####Choose the data####

dataset=(S1toS2 S1toS2_windows)

data_name=(S1toS2_exact_match_55_20kb.bed S1toS2_exact_match_100kb_windowed.bed)

combine_options=("-score T -strand T" "-strand T")

##### 