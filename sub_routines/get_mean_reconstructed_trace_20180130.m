function [x_rec,rec_trace] = get_mean_reconstructed_trace_20180130(...
    bool_3, bool_5, x, x_3, x_5, isoform_levels)
% generates the reconstructed mean traces from the isoform abundance for
% display (overlayed on the Rend-seq data).
%
% Input:
% bool_3: bool_3(i) = 1 if x(i) is a 3' end, 0 otherwise.
% bool_5: bool_5(i) = 1 if x(i) is a 5' end, 0 otherwise.
% x: positions of ends in roi.
% x_3: positions of 3' ends in roi.
% x_5: positions of 5' ends in roi.
% isoform_levels: isoform level array, entry at (i,j) equals abundance of
%                 mRNA isoform starting at x_5(i) and ending at x_3(j).
%
% Outputs:
% x_rec: positions for reconstructed trace.
% rec_trace: reconstructed trace (mean levels between ends). 


% initialization
n_ends = length(x);
x_rec = x(1):x(end);
rec_trace = zeros(max(x),1);

% looping through regions, summing all isoforms overlapping with region
for i = 2:n_ends
    if bool_5(i)
        ind3 = find(x_3>x(i));
        ind5 = find(x_5<x(i));
        rec_level = sum(sum(isoform_levels(ind5,ind3)));
    elseif bool_3(i)
        ind3 = find(x(i)==x_3);
        ind5 = find(x_5<x(i));
        rec_level = sum(sum(isoform_levels(ind5,ind3:end)));
    end
    
    % allocating 
    rec_trace(x(i-1):x(i)) = rec_level;
    
end

% trimming unnecesary region.
ind_range = 1:length(rec_trace);
rec_trace = rec_trace((ind_range>=min(x))' & (ind_range<=max(x))');
