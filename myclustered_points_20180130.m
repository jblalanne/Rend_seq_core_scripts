function [clustered_points_ind, clustered_points] = myclustered_points_20180130(points,distance_threshold)
% returns clusters of points (1D positions) within a given distance
% threshold.
%
% Inputs: 
% points: list of positions (1D).
% distance_threshold: distance for clustering
%
% Outputs:
% clustered_points_id: data structure with indices of clustered points
% clustered_points: data structure with clustered points


% initialization
clustered_points_ind = [];
counter_cluster = 1;
new_cluster = 1;
n_points = length(points);

% boolean for distance between points
close_points_bool = diff(points)<distance_threshold;

% going through list and grouping points based on distance
for i = 1:n_points
    
    if new_cluster
        clustered_points_ind{counter_cluster} = i;
        new_cluster = 0;
    else
        clustered_points_ind{counter_cluster} = [clustered_points_ind{counter_cluster}; i];
    end
    
    % moving to a next cluster if farther than the distance threshold.
    if i < n_points
        if ~close_points_bool(i)
            counter_cluster = counter_cluster+1;
            new_cluster = 1;
        end
    end
end

for i = 1:length(clustered_points_ind)
    clustered_points{i} = points(clustered_points_ind{i});
end














