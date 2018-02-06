function [isoform_levels, x_rec, rec_trace] = ...
    isoform_quantification_20180130(x_5,x_3,data_3)
% Computes the mRNA isoforms from read count data and ends positions, and
% returns the reconstructed mean trace for display.
%
% Input:
% x_5: positions of 5' ends identified in roi.
% x_3: positions of 3' ends identified in roi.
% data_3: counts for 3' ends of mapped reads in roi.
%
% Outputs: 
% isoform_levels: isoform level array, entry at (i,j) equals abundance of
%                 mRNA isoform starting at x_5(i) and ending at x_3(j).
%                 Units is in reads/nt.
% x_rec: positions for reconstructed trace.
% rec_trace: reconstructed trace (mean levels between ends). 


%% sorting the ends
x_3 = sort(x_3,'ascend');
x_5 = sort(x_5,'ascend');

%% re-format of positions (ordered irrespective of type)
[x,bool_3,bool_5] = reorder_positions_20180130(x_3,x_5);
    
%% quantify levels between ends
gap_small = 3;
gap_read_small = 14;
gap_read_large = 45;
[levels_subregions,lengths_subregions] = ...
    get_level_between_ends_20180130(data_3, x, bool_3, bool_5,...
    gap_small, gap_read_small, gap_read_large);

%% get f and rho 
% iteratively get f and rho (effectively 3' ends with negative apparent
% readthrough, or 5' ends with negative promoter strength). 
% See Methods in publication for description of f and rho and derivation.
max_iter = 4;
[f, rho] = get_f_and_rho_iterative_20180130(...
    x, x_5, bool_3, bool_5, levels_subregions, lengths_subregions, max_iter);

%% get isoform levels
isoform_levels = get_isoform_levels_20180130(x,x_5,x_3,f,rho);

%% mean reconstructed trace
[x_rec,rec_trace] = get_mean_reconstructed_trace_20180130(...
    bool_3, bool_5, x, x_3, x_5, isoform_levels);
