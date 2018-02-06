function quantify_mRNA_isoforms_20180130(data)



%% get ends in region
get_ends_roi_20180130();

%% get mRNA isoform levels
[isoform_levels, x_rec, rec_trace] = ...
    isoform_quantification_20180130(x_5,x_3,data_3);


[levels,termination_probabilities,promoter_strength,...
    isoform_levels, x_sim, simulated_trace,data_3pr,data_5pr, x_3, x_5] =...
    get_isoforms_v3(gap_small,gap_read_small,gap_read_large,buffer_operon,...
    z_thresh,width_z,z,data,start_genes,end_genes,strand_operon,...
    known_5pr,known_3pr,x5_rem,x3_rem)
