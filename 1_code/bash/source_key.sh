#!/bin/bash

# Switch on/off

Key_pre_analysis=false # from "directory for files" to "bedGraph to BigWig"
if ${Key_pre_analysis}; then
    Key_cutadapt=true
    Key_bowtie2=true
    Key_sort_count=true
    Key_sep_plus_minus=true # contains also calculations of length distribution and dinucleotide content tables.
    Key_bedgraph_BigWig=true
    Key_TS_NTS=true
    
fi

Key_downstream_analysis=true # from "aligning to datasets" to "combining both strands in a file"
if ${Key_downstream_analysis}; then
    Key_alignment=false
    Key_combineWindows=true
    Key_moreInfo=true
    Key_rpkm=true
    Key_file=true

fi

