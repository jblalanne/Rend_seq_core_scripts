function [x_3, x_5, isoform_levels] = quantify_plot_isoforms_20180130(start_region,stop_region,...
    strand_region,z_thresh,width_z,z_peak_oi,data_oi,known_5,known_3,...
    spurious_5,spurious_3,genome_size,annotation_file,annotation_dir)
% identifies and quantifies mRNA isoforms from Rend-seq data for region of 
% interest. Also displays the data and isoforms.
%
% Inputs: 
% start_region: start of region of interest (position in nt on chromosome).
%           On the reverse strand, start is the most downstream position
%           (i.e., start_region < stop_region even on reverse strand).
% stop_region: end of region of interest (position in nt on chromosome).
%           (stop_region > start_region even on reverse strand). 
% strand_region: strand of region of interest (1: forward, 0 reverse). 
% z_thresh: threshold on peak z score for end identification.
% width_z: distance (nt) within which multiple identified ends are merged. 
% z_peak_oi: peak z score for dataset of interest (genome size by 4 array).
% data_oi: Rend-seq read counts for dataset of interest (genome size by 4 array).
% known_5: known 5' ends missed by automated analysis (chromosomal
%           coordinates).
% known_3: known 3' ends missed by automated analysis (chromosomal
%           coordinates).
% spurious_5: spuriously identified 5' ends to be removed (chromosomal
%           coordinates).
% spurious_3: spuriously identified 3' ends to be removed (chromosomal
%           coordinates).
% genome_size: size of genome (in nt).
% annotation_file: name of GenBank annotation file. 
% annotation_dir: full path of directory containing annotation file.
%
% Outputs:
% x_3: position of 3' ends in region of interest.
% x_5: position of 5' ends in region of interest.
% isoform_levels: isoform level array, entry at (i,j) equals abundance of
%                 mRNA isoform starting at x_5(i) and ending at x_3(j).
%                 Units is in reads/nt.

%% get ends in region
[data_3,data_5, x_3, x_5] =...
    get_ends_roi_20180130(start_region,stop_region,strand_region,...
    z_thresh,width_z,z_peak_oi,data_oi,known_5,known_3,spurious_5,spurious_3);


%% get mRNA isoform levels
[isoform_levels, x_rec, rec_trace] = ...
    isoform_quantification_20180130(x_5,x_3,data_3);


%% plotting
% Rend-seq traces with reconstructed traces and isoforms with relative
% abundance.
plot_isoforms_20180131(x_rec,rec_trace,x_5,x_3,data_5,data_3,strand_region,...
    isoform_levels,genome_size,annotation_file,annotation_dir)

% convert back x_3 and x_5 to chromosome coordinate if on the reverse
% strand
if ~strand_region
    x_5 = genome_size-x_5+1;
    x_3 = genome_size-x_3+1;
end


