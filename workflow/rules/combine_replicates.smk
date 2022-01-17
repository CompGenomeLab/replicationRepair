rule combine_replicates:
    input:
        HXA64A_plus="resources/samples/XR/HXA64A1_ATCACG_hg19_sorted_plus.bed",
        HXA64A_minus="resources/samples/XR/HXA64A1_ATCACG_hg19_sorted_minus.bed",
        HXA64B_plus="resources/samples/XR/HXA64B7_CAGATC_hg19_sorted_plus.bed",
        HXA64B_minus="resources/samples/XR/HXA64B7_CAGATC_hg19_sorted_minus.bed",
        HXACA_plus="resources/samples/XR/HXACA4_TGACCA_hg19_sorted_plus.bed",
        HXACA_minus="resources/samples/XR/HXACA4_TGACCA_hg19_sorted_minus.bed",
        HXACB_plus="resources/samples/XR/HXACB10_TAGCTT_hg19_sorted_plus.bed",
        HXACB_minus="resources/samples/XR/HXACB10_TAGCTT_hg19_sorted_minus.bed",
        HXE64A_plus="resources/samples/XR/HXE64A2_CGATGT_hg19_sorted_plus.bed",
        HXE64A_minus="resources/samples/XR/HXE64A2_CGATGT_hg19_sorted_minus.bed",
        HXE64B_plus="resources/samples/XR/HXE64B8_ACTTGA_hg19_sorted_plus.bed",
        HXE64B_minus="resources/samples/XR/HXE64B8_ACTTGA_hg19_sorted_minus.bed", 
        HXECA_plus="resources/samples/XR/HXECA5_ACAGTG_hg19_sorted_plus.bed",
        HXECA_minus="resources/samples/XR/HXECA5_ACAGTG_hg19_sorted_minus.bed",
        HXECB_plus="resources/samples/XR/HXECB11_GGCTAC_hg19_sorted_plus.bed",
        HXECB_minus="resources/samples/XR/HXECB11_GGCTAC_hg19_sorted_minus.bed",        
        HXL64A_plus="resources/samples/XR/HXL64A3_TTAGGC_hg19_sorted_plus.bed",
        HXL64A_minus="resources/samples/XR/HXL64A3_TTAGGC_hg19_sorted_minus.bed",
        HXL64B_plus="resources/samples/XR/HXL64B9_GATCAG_hg19_sorted_plus.bed",
        HXL64B_minus="resources/samples/XR/HXL64B9_GATCAG_hg19_sorted_minus.bed",   
        HXLCA_plus="resources/samples/XR/HXLCA6_GCCAAT_hg19_sorted_plus.bed",
        HXLCA_minus="resources/samples/XR/HXLCA6_GCCAAT_hg19_sorted_minus.bed",
        HXLCB_plus="resources/samples/XR/HXLCB12_CTTGTA_hg19_sorted_plus.bed",
        HXLCB_minus="resources/samples/XR/HXLCB12_CTTGTA_hg19_sorted_minus.bed",      
        Hela15X2_plus="resources/samples/XR/R19026421-2019XR1-Hela15X2_combined_R1_hg19_sorted_plus.bed",
        Hela15X2_minus="resources/samples/XR/R19026421-2019XR1-Hela15X2_combined_R1_hg19_sorted_minus.bed",
        Hela15X7_plus="resources/samples/XR/R19033030-2019XR3-Hela15X7_combined_R1_hg19_sorted_plus.bed",
        Hela15X7_minus="resources/samples/XR/R19033030-2019XR3-Hela15X7_combined_R1_hg19_sorted_minus.bed",
        Hela35X3_plus="resources/samples/XR/R19026421-2019XR1-Hela35X3_combined_R1_hg19_sorted_plus.bed",
        Hela35X3_minus="resources/samples/XR/R19026421-2019XR1-Hela35X3_combined_R1_hg19_sorted_minus.bed",   
        Hela35X8_plus="resources/samples/XR/R19033030-2019XR3-Hela35X8_combined_R1_hg19_sorted_plus.bed",
        Hela35X8_minus="resources/samples/XR/R19033030-2019XR3-Hela35X8_combined_R1_hg19_sorted_minus.bed", 
        HDA64A_plus="resources/samples/DS/HDA64A1_ATCACG_hg19_sorted_ds_dipyrimidines_plus.bed",
        HDA64A_minus="resources/samples/DS/HDA64A1_ATCACG_hg19_sorted_ds_dipyrimidines_minus.bed",
        HDA64B_plus="resources/samples/DS/HDA64B19_GTGAAA_hg19_sorted_ds_dipyrimidines_plus.bed",
        HDA64B_minus="resources/samples/DS/HDA64B19_GTGAAA_hg19_sorted_ds_dipyrimidines_minus.bed",
        HDACA_plus="resources/samples/DS/HDACA6_GCCAAT_hg19_sorted_ds_dipyrimidines_plus.bed",
        HDACA_minus="resources/samples/DS/HDACA6_GCCAAT_hg19_sorted_ds_dipyrimidines_minus.bed",
        HDACB_plus="resources/samples/DS/HDACB23_GAGTGG_hg19_sorted_ds_dipyrimidines_plus.bed",
        HDACB_minus="resources/samples/DS/HDACB23_GAGTGG_hg19_sorted_ds_dipyrimidines_minus.bed",
        HDE64A_plus="resources/samples/DS/HDE64A4_TGACCA_hg19_sorted_ds_dipyrimidines_plus.bed",
        HDE64A_minus="resources/samples/DS/HDE64A4_TGACCA_hg19_sorted_ds_dipyrimidines_minus.bed",
        HDE64B_plus="resources/samples/DS/HDE64B20_GTGGCC_hg19_sorted_ds_dipyrimidines_plus.bed",
        HDE64B_minus="resources/samples/DS/HDE64B20_GTGGCC_hg19_sorted_ds_dipyrimidines_minus.bed", 
        HDECA_plus="resources/samples/DS/HDECA10_TAGCTT_hg19_sorted_ds_dipyrimidines_plus.bed",
        HDECA_minus="resources/samples/DS/HDECA10_TAGCTT_hg19_sorted_ds_dipyrimidines_minus.bed",
        HDECB_plus="resources/samples/DS/HDECB25_ACTGAT_hg19_sorted_ds_dipyrimidines_plus.bed",
        HDECB_minus="resources/samples/DS/HDECB25_ACTGAT_hg19_sorted_ds_dipyrimidines_minus.bed",        
        HDL64A_plus="resources/samples/DS/HDL64A5_ACAGTG_hg19_sorted_ds_dipyrimidines_plus.bed",
        HDL64A_minus="resources/samples/DS/HDL64A5_ACAGTG_hg19_sorted_ds_dipyrimidines_minus.bed",
        HDL64B_plus="resources/samples/DS/HDL64B22_CGTACG_hg19_sorted_ds_dipyrimidines_plus.bed",
        HDL64B_minus="resources/samples/DS/HDL64B22_CGTACG_hg19_sorted_ds_dipyrimidines_minus.bed",   
        HDLCA_plus="resources/samples/DS/HDLCA12_CTTGTA_hg19_sorted_ds_dipyrimidines_plus.bed",
        HDLCA_minus="resources/samples/DS/HDLCA12_CTTGTA_hg19_sorted_ds_dipyrimidines_minus.bed",
        HDLCB_plus="resources/samples/DS/HDLCB27_ATTCCT_hg19_sorted_ds_dipyrimidines_plus.bed",
        HDLCB_minus="resources/samples/DS/HDLCB27_ATTCCT_hg19_sorted_ds_dipyrimidines_minus.bed",      
        CT5H2_HelaD1_plus="resources/samples/DS/NT4H1-CT5H2-HelaD1-5R2h2_combined_hg19_sorted_ds_dipyrimidines_plus.bed",
        CT5H2_HelaD1_minus="resources/samples/DS/NT4H1-CT5H2-HelaD1-5R2h2_combined_R1_hg19_sorted_ds_dipyrimidines_minus.bed",
        HD2TCD4_HelaD1_plus="resources/samples/DS/R19029847-HD2TCD4-HelaD1-5R2h1_combined_hg19_sorted_ds_dipyrimidines_plus.bed",
        HD2TCD4_HelaD1_minus="resources/samples/DS/R19029847-HD2TCD4-HelaD1-5R2h1_combined_hg19_sorted_ds_dipyrimidines_minus.bed",
        CT5H2_HelaD3_plus="resources/samples/DS/XT4H1-CT5H2-HelaD3-5R2h2_combined_hg19_sorted_ds_dipyrimidines_plus.bed",
        CT5H2_HelaD3_minus="resources/samples/DS/XT4H1-CT5H2-HelaD3-5R2h2_combined_hg19_sorted_ds_dipyrimidines_minus.bed",   
        HD2TCD4_HelaD3_plus="resources/samples/DS/R19029847-HD2TCD4-HelaD3-5R2h1_combined_hg19_sorted_ds_dipyrimidines_plus.bed",
        HD2TCD4_HelaD3_minus="resources/samples/DS/R19029847-HD2TCD4-HelaD3-5R2h1_combined_hg19_sorted_ds_dipyrimidines_minus.bed", 
    output:
        HXA64_minus="resources/samples/XR/HXA64_comb_hg19_sorted_minus.bed",
        HXA64_plus="resources/samples/XR/HXA64_comb_hg19_sorted_plus.bed",
        HXAC_minus="resources/samples/XR/HXAC_comb_hg19_sorted_minus.bed",
        HXAC_plus="resources/samples/XR/HXAC_comb_hg19_sorted_plus.bed",
        HXE64_minus="resources/samples/XR/HXE64_comb_hg19_sorted_minus.bed",
        HXE64_plus="resources/samples/XR/HXE64_comb_hg19_sorted_plus.bed",
        HXEC_minus="resources/samples/XR/HXEC_comb_hg19_sorted_minus.bed",
        HXEC_plus="resources/samples/XR/HXEC_comb_hg19_sorted_plus.bed",
        HXL64_minus="resources/samples/XR/HXL64_comb_hg19_sorted_minus.bed",
        HXL64_plus="resources/samples/XR/HXL64_comb_hg19_sorted_plus.bed",
        HXLC_minus="resources/samples/XR/HXLC_comb_hg19_sorted_minus.bed",
        HXLC_plus="resources/samples/XR/HXLC_comb_hg19_sorted_plus.bed",
        Hela15X2_plus="resources/samples/XR/Hela15X2_comb_hg19_sorted_plus.bed",
        Hela15X2_minus="resources/samples/XR/Hela15X2_comb_hg19_sorted_minus.bed",
        Hela35X3_plus="resources/samples/XR/Hela35X3_comb_hg19_sorted_plus.bed",
        Hela35X3_minus="resources/samples/XR/Hela35X3_comb_hg19_sorted_minus.bed",
        HDA64_minus="resources/samples/DS/HDA64_comb_hg19_sorted_ds_dipyrimidines_minus.bed",
        HDA64_plus="resources/samples/DS/HDA64_comb_hg19_sorted_ds_dipyrimidines_plus.bed",
        HDAC_minus="resources/samples/DS/HDAC_comb_hg19_sorted_ds_dipyrimidines_minus.bed",
        HDAC_plus="resources/samples/DS/HDAC_comb_hg19_sorted_ds_dipyrimidines_plus.bed",
        HDE64_minus="resources/samples/DS/HDE64_comb_hg19_sorted_ds_dipyrimidines_minus.bed",
        HDE64_plus="resources/samples/DS/HDE64_comb_hg19_sorted_ds_dipyrimidines_plus.bed",
        HDEC_minus="resources/samples/DS/HDEC_comb_hg19_sorted_ds_dipyrimidines_minus.bed",
        HDEC_plus="resources/samples/DS/HDEC_comb_hg19_sorted_ds_dipyrimidines_plus.bed",
        HDL64_minus="resources/samples/DS/HDL64_comb_hg19_sorted_ds_dipyrimidines_minus.bed",
        HDL64_plus="resources/samples/DS/HDL64_comb_hg19_sorted_ds_dipyrimidines_plus.bed",
        HDLC_minus="resources/samples/DS/HDLC_comb_hg19_sorted_ds_dipyrimidines_minus.bed",
        HDLC_plus="resources/samples/DS/HDLC_comb_hg19_sorted_ds_dipyrimidines_plus.bed",
        HelaD1_5R2h2_plus="resources/samples/DS/HelaD1_5R2h2_comb_hg19_sorted_ds_dipyrimidines_plus.bed",
        HelaD1_5R2h2_minus="resources/samples/DS/HelaD1_5R2h2_comb_hg19_sorted_ds_dipyrimidines_minus.bed",
        HelaD3_5R2h2_plus="resources/samples/DS/HelaD3_5R2h2_comb_hg19_sorted_ds_dipyrimidines_plus.bed",
        HelaD3_5R2h2_minus="resources/samples/DS/HelaD3_5R2h2_comb_hg19_sorted_ds_dipyrimidines_minus.bed",
    log:
        "logs/combine_replicates.log",
    benchmark:
        "logs/combine_replicates.benchmark.txt",
    shell:  
        """
        (echo "`date -R`: Combining all replicates..." &&
        cat {input.HXA64A_plus} {input.HXA64B_plus} > {output.HXA64_plus} &&
        cat {input.HXA64A_minus} {input.HXA64B_minus} > {output.HXA64_minus} &&
        cat {input.HXACA_plus} {input.HXACB_plus} > {output.HXAC_plus} &&
        cat {input.HXACA_minus} {input.HXACB_minus} > {output.HXAC_minus} &&
        cat {input.HXE64A_plus} {input.HXE64B_plus} > {output.HXE64_plus} &&
        cat {input.HXE64A_minus} {input.HXE64B_minus} > {output.HXE64_minus} &&
        cat {input.HXECA_plus} {input.HXECB_plus} > {output.HXEC_plus} &&
        cat {input.HXECA_minus} {input.HXECB_minus} > {output.HXEC_minus} &&
        cat {input.HXL64A_plus} {input.HXL64B_plus} > {output.HXL64_plus} &&
        cat {input.HXL64A_minus} {input.HXL64B_minus} > {output.HXL64_minus} &&
        cat {input.HXLCA_plus} {input.HXLCB_plus} > {output.HXLC_plus} &&
        cat {input.HXLCA_minus} {input.HXLCB_minus} > {output.HXLC_minus} &&
        cat {input.Hela15X2_plus} {input.Hela15X7_plus} > {output.Hela15X2_plus} &&
        cat {input.Hela15X2_minus} {input.Hela15X7_minus} > {output.Hela15X2_minus} &&
        cat {input.Hela35X3_plus} {input.Hela35X8_plus} > {output.Hela35X3_plus} &&
        cat {input.Hela35X3_minus} {input.Hela35X8_minus} > {output.Hela35X3_minus} &&
        cat {input.HDA64A_plus} {input.HDA64B_plus} > {output.HDA64_plus} &&
        cat {input.HDA64A_minus} {input.HDA64B_minus} > {output.HDA64_minus} &&
        cat {input.HDACA_plus} {input.HDACB_plus} > {output.HDAC_plus} &&
        cat {input.HDACA_minus} {input.HDACB_minus} > {output.HDAC_minus} &&
        cat {input.HDE64A_plus} {input.HDE64B_plus} > {output.HDE64_plus} &&
        cat {input.HDE64A_minus} {input.HDE64B_minus} > {output.HDE64_minus} &&
        cat {input.HDECA_plus} {input.HDECB_plus} > {output.HDEC_plus} &&
        cat {input.HDECA_minus} {input.HDECB_minus} > {output.HDEC_minus} &&
        cat {input.HDL64A_plus} {input.HDL64B_plus} > {output.HDL64_plus} &&
        cat {input.HDL64A_minus} {input.HDL64B_minus} > {output.HDL64_minus} &&
        cat {input.HDLCA_plus} {input.HDLCB_plus} > {output.HDLC_plus} &&
        cat {input.HDLCA_minus} {input.HDLCB_minus} > {output.HDLC_minus} &&
        cat {input.CT5H2_HelaD1_plus} {input.HD2TCD4_HelaD1_plus} > {output.HelaD1_5R2h2_plus} &&
        cat {input.CT5H2_HelaD1_minus} {input.HD2TCD4_HelaD1_minus} > {output.HelaD1_5R2h2_minus} &&
        cat {input.CT5H2_HelaD3_plus} {input.HD2TCD4_HelaD3_plus} > {output.HelaD3_5R2h2_plus} &&
        cat {input.CT5H2_HelaD3_minus} {input.HD2TCD4_HelaD3_minus} > {output.HelaD3_5R2h2_minus} &&
        echo "`date -R`: Success! All replicates are combined." || 
        {{ echo "`date -R`: Process failed..."; exit 1; }} ) > {log} 2>&1
        """

rule combine_replicates_sim:
    input:
        HXA64A="resources/samples/sim/HXA64A1_ATCACG_hg19_xr_sim.bed",
        HXA64B="resources/samples/sim/HXA64B7_CAGATC_hg19_xr_sim.bed",
        HXACA="resources/samples/sim/HXACA4_TGACCA_hg19_xr_sim.bed",
        HXACB="resources/samples/sim/HXACB10_TAGCTT_hg19_xr_sim.bed",
        HXE64A="resources/samples/sim/HXE64A2_CGATGT_hg19_xr_sim.bed",
        HXE64B="resources/samples/sim/HXE64B8_ACTTGA_hg19_xr_sim.bed",
        HXECA="resources/samples/sim/HXECA5_ACAGTG_hg19_xr_sim.bed",
        HXECB="resources/samples/sim/HXECB11_GGCTAC_hg19_xr_sim.bed",        
        HXL64A="resources/samples/sim/HXL64A3_TTAGGC_hg19_xr_sim.bed",
        HXL64B="resources/samples/sim/HXL64B9_GATCAG_hg19_xr_sim.bed",  
        HXLCA="resources/samples/sim/HXLCA6_GCCAAT_hg19_xr_sim.bed",
        HXLCB="resources/samples/sim/HXLCB12_CTTGTA_hg19_xr_sim.bed",     
        Hela15X2="resources/samples/sim/R19026421-2019XR1-Hela15X2_combined_R1_hg19_xr_sim.bed",
        Hela15X7="resources/samples/sim/R19033030-2019XR3-Hela15X7_combined_R1_hg19_xr_sim.bed",
        Hela35X3="resources/samples/sim/R19026421-2019XR1-Hela35X3_combined_R1_hg19_xr_sim.bed",  
        Hela35X8="resources/samples/sim/R19033030-2019XR3-Hela35X8_combined_R1_hg19_xr_sim.bed",
        HDA64A="resources/samples/sim/HDA64A1_ATCACG_hg19_ds_sim.bed",
        HDA64B="resources/samples/sim/HDA64B19_GTGAAA_hg19_ds_sim.bed",
        HDACA="resources/samples/sim/HDACA6_GCCAAT_hg19_ds_sim.bed",
        HDACB="resources/samples/sim/HDACB23_GAGTGG_hg19_ds_sim.bed",
        HDE64A="resources/samples/sim/HDE64A4_TGACCA_hg19_ds_sim.bed",
        HDE64B="resources/samples/sim/HDE64B20_GTGGCC_hg19_ds_sim.bed",
        HDECA="resources/samples/sim/HDECA10_TAGCTT_hg19_ds_sim.bed",
        HDECB="resources/samples/sim/HDECB25_ACTGAT_hg19_ds_sim.bed",        
        HDL64A="resources/samples/sim/HDL64A5_ACAGTG_hg19_ds_sim.bed",
        HDL64B="resources/samples/sim/HDL64B22_CGTACG_hg19_ds_sim.bed",  
        HDLCA="resources/samples/sim/HDLCA12_CTTGTA_hg19_ds_sim.bed",
        HDLCB="resources/samples/sim/HDLCB27_ATTCCT_hg19_ds_sim.bed",     
        HD2TCD4_HelaD1="resources/samples/sim/R19029847-HD2TCD4-HelaD1-5R2h1_combined_hg19_ds_sim.bed",
        CT5H2_HelaD1="resources/samples/sim/NT4H1-CT5H2-HelaD1-5R2h2_combined_hg19_ds_sim.bed",
        HD2TCD4_HelaD3="resources/samples/sim/R19029847-HD2TCD4-HelaD3-5R2h1_combined_hg19_ds_sim.bed",  
        CT5H2_HelaD3="resources/samples/sim/XT4H1-CT5H2-HelaD3-5R2h2_combined_hg19_ds_sim.bed",
    output:
        HXA64="resources/samples/sim/HXA64_comb_hg19_xr_sim.bed",
        HXAC="resources/samples/sim/HXAC_comb_hg19_xr_sim.bed",
        HXE64="resources/samples/sim/HXE64_comb_hg19_xr_sim.bed",
        HXEC="resources/samples/sim/HXEC_comb_hg19_xr_sim.bed",
        HXL64="resources/samples/sim/HXL64_comb_hg19_xr_sim.bed",
        HXLC="resources/samples/sim/HXLC_comb_hg19_xr_sim.bed",
        Hela15X2="resources/samples/sim/Hela15X2_comb_hg19_xr_sim.bed",
        Hela35X3="resources/samples/sim/Hela35X3_comb_hg19_xr_sim.bed",
        HDA64="resources/samples/sim/HDA64_comb_hg19_ds_sim.bed",
        HDAC="resources/samples/sim/HDAC_comb_hg19_ds_sim.bed",
        HDE64="resources/samples/sim/HDE64_comb_hg19_ds_sim.bed",
        HDEC="resources/samples/sim/HDEC_comb_hg19_ds_sim.bed",      
        HDL64="resources/samples/sim/HDL64_comb_hg19_ds_sim.bed",  
        HDLC="resources/samples/sim/HDLC_comb_hg19_ds_sim.bed",    
        HelaD1_5R2h2="resources/samples/sim/HelaD1_5R2h2_comb_hg19_ds_sim.bed",
        HelaD3_5R2h2="resources/samples/sim/HelaD3_5R2h2_comb_hg19_ds_sim.bed",
    log:
        "logs/combine_replicates_sim.log",
    benchmark:
        "logs/combine_replicates_sim.benchmark.txt",
    shell:  
        """
        (echo "`date -R`: Combining all replicates..." &&
        cat {input.HXA64A} {input.HXA64B} > {output.HXA64} &&
        cat {input.HXACA} {input.HXACB} > {output.HXAC} &&
        cat {input.HXE64A} {input.HXE64B} > {output.HXE64} &&
        cat {input.HXECA} {input.HXECB} > {output.HXEC} &&
        cat {input.HXL64A} {input.HXL64B} > {output.HXL64} &&
        cat {input.HXLCA} {input.HXLCB} > {output.HXLC} &&
        cat {input.Hela15X2} {input.Hela15X7} > {output.Hela15X2} &&
        cat {input.Hela35X3} {input.Hela35X8} > {output.Hela35X3} &&
        cat {input.HDA64A} {input.HDA64B} > {output.HDA64} &&
        cat {input.HDACA} {input.HDACB} > {output.HDAC} &&
        cat {input.HDE64A} {input.HDE64B} > {output.HDE64} &&
        cat {input.HDECA} {input.HDECB} > {output.HDEC} &&
        cat {input.HDL64A} {input.HDL64B} > {output.HDL64} &&
        cat {input.HDLCA} {input.HDLCB} > {output.HDLC} &&
        cat {input.HD2TCD4_HelaD1} {input.CT5H2_HelaD1} > {output.HelaD1_5R2h2} &&
        cat {input.HD2TCD4_HelaD3} {input.CT5H2_HelaD3} > {output.HelaD3_5R2h2} &&
        echo "`date -R`: Success! All replicates are combined." || 
        {{ echo "`date -R`: Process failed..."; exit 1; }} ) > {log} 2>&1
        """

rule combine_methyl:
    input:
        CpG="results/methyl/GSM3633947_ENCFF696OLO_methylation_state_at_CpG/GSM3633947_ENCFF696OLO_methylation_state_at_CpG_hg19_org_intergenic.bed",
        CHH="results/methyl/GSM3633947_ENCFF474VQX_methylation_state_at_CHH/GSM3633947_ENCFF474VQX_methylation_state_at_CHH_hg19_org_intergenic.bed",
        CHG="results/methyl/GSM3633947_ENCFF093WZD_methylation_state_at_CHG/GSM3633947_ENCFF093WZD_methylation_state_at_CHG_hg19_org_intergenic.bed",
    output:
        CpG="results/methyl/GSM3633947_ENCFF696OLO_methylation_state_at_CpG/GSM3633947_ENCFF696OLO_methylation_state_at_CpG_hg19_org_intergenic_shuf1m.bed",
        CHH="results/methyl/GSM3633947_ENCFF474VQX_methylation_state_at_CHH/GSM3633947_ENCFF474VQX_methylation_state_at_CHH_hg19_org_intergenic_shuf1m.bed",
        CHG="results/methyl/GSM3633947_ENCFF093WZD_methylation_state_at_CHG/GSM3633947_ENCFF093WZD_methylation_state_at_CHG_hg19_org_intergenic_shuf1m.bed",
        comb="results/regions/methylation_shuf_1m_windows_1_3.bed",
    log:
        "logs/combine_methyl.log",
    benchmark:
        "logs/combine_methyl.benchmark.txt",
    shell:  
        """
        (echo "`date -R`: Getting random million reads for each file..." &&
        shuf -n 1000000 {input.CpG} | awk '{{print $1"\\t"$2"\\t"$3"\\t""CpG_1""\\t"$5"\\t"$6}}' > {output.CpG} &&
        shuf -n 1000000 {input.CHG} | awk '{{print $1"\\t"$2"\\t"$3"\\t""CHG_1""\\t"$5"\\t"$6}}' > {output.CHG} &&
        shuf -n 1000000 {input.CHH} | awk '{{print $1"\\t"$2"\\t"$3"\\t""CHH_1""\\t"$5"\\t"$6}}' > {output.CHH} &&
        echo "`date -R`: Success! All files are combined." || 
        {{ echo "`date -R`: Process failed..."; exit 1; }} ) > {log} 2>&1
        
        (echo "`date -R`: Combining methylation files..." &&
        cat {output.CpG} {output.CHH} {output.CHG}  > {output.comb} &&
        echo "`date -R`: Success! All files are combined." || 
        {{ echo "`date -R`: Process failed..."; exit 1; }} ) >> {log} 2>&1
        """