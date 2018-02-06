function peak_positions = get_z_peaks_20180130(z_peak,z_thresh,d_thresh)
% returns positions with z score above threshold (returns maximum z score
% positions in d_thresh neighborhood if multiple close by).
%
% Inputs:
% z_peak: peak z score vs. position. 3D array:
%           dimension 1: different datasets.
%           dimension 2: genome position
%           dimension 3: type of reads: 3f, 3r, 5f, 5r.
% z_thresh: threshold for the peak z score (peaks have peak z > z_thresh)
% d_thresh: distance threshold for grouping neighboring peaks.
%
% Outputs:
% peak_positions: peak positions for:
%           dimension 1: different datasets
%           dimension 2: type of ends (3'f, 3'r, 5'f, 5'r)
%


% get all positions above z score threshold
all_peak_positions = [];
for i = 1:size(z_peak,1)
    for j = 1:4
        all_peak_positions{i,j} = find((z_peak(i,:,j)>z_thresh));
    end
end

% find peaks that are closer than dist_thresh, group them.
clustered_points = [];
for i = 1:size(z_peak,1)
    for j = 1:4
        [~, clustered_points{i,j}] = myclustered_points_20180130(all_peak_positions{i,j},d_thresh);
    end
end

% selecting the peak within each cluster of position with the largest z score. 
% takes the first position if two have identical z score (unlikely).
peak_positions = [];
for i = 1:size(z_peak,1)
    for j = 1:4
        type_specific_peaks = [];
        for k = 1:length(clustered_points{i,j})
            cluster = clustered_points{i,j}{k};
            ind = find(z_peak(i,cluster,j)==...
                max(z_peak(i,cluster,j)),1,'first');
            
            type_specific_peaks = [type_specific_peaks cluster(ind)];
        end
        peak_positions{i,j} = type_specific_peaks;
    end
end