
sminus <- filter(pB1_data, repdomains == "DTZ", sample_strand == "-", windows > 0)
sminus <- sminus$xr_ds
splus <- filter(pB1_data, repdomains == "DTZ", sample_strand == "+", windows > 0)
splus <- splus$xr_ds
ks.test(sminus, splus)
