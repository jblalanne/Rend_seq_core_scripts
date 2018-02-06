function [levels,termination_probabilities,promoter_strength,...
    isoform_levels, x_sim, simulated_trace,data_3pr,data_5pr, x_3, x_5] =...
    get_isoforms_20180130(gap_small,gap_read_small,gap_read_large,buffer,...
    z_thresh,width_z,z_peak,data,starts_genes,stops_genes,strand_region,...
    known_5pr,known_3pr,x5_rem,x3_rem)


if strand_region
    start_range = min([known_5pr min(starts_genes)])-buffer;
    end_range = max([known_3pr max(starts_genes)])+buffer;
else
    start_range = min([known_3pr min(starts_genes)])-buffer;
    end_range = max([known_5pr max(starts_genes)])+buffer;
end
range = start_range:end_range;


% looking for 3' and 5' peaks in roi
if strand_region(1)
    ends_3pr = find(z_peak(1,:,1)>z_thresh);
    ends_5pr = find(z_peak(1,:,3)>z_thresh);
else
    ends_3pr = find(z_peak(1,:,2)>z_thresh);
    ends_5pr = find(z_peak(1,:,4)>z_thresh);    
end

% keep only the ends in roi
ends_3pr = ends_3pr(ends_3pr>min(range) & ends_3pr<max(range));
ends_5pr = ends_5pr(ends_5pr>min(range) & ends_5pr<max(range));

% remove redundant positions
temp = ends_3pr;
temp(diff(temp)<=width_z) = [];
ends_3pr = temp;
temp = ends_5pr;
temp(diff(temp)<=width_z) = [];
ends_5pr = temp;

% adding known ends
ends_5pr = sort([ends_5pr known_5pr]);
ends_3pr = sort([ends_3pr known_3pr]);


% get rid of 5' ends after the last gene and the 3' ends before the first
% gene
if strand_region
    
    if isempty(ends_5pr)
        ends_5pr(end+1) = min(starts_genes);
    elseif min(ends_5pr)>min(starts_genes)
        % no 5' end before the start of the gene in the buffer region.
        % add an "artificial" 5' end at the beginning of the region of
        % interest.
        ends_5pr(end+1) = min(starts_genes);
    end
    
    if isempty(ends_3pr)
        ends_3pr(end+1) = max(stops_genes);
    elseif max(ends_3pr)<max(stops_genes)
        ends_3pr(end+1) = max(stops_genes);
    end
    
    ends_3pr(ends_3pr<min(starts_genes)) = [];
    ends_5pr(ends_5pr>max(stops_genes)) = [];
    
else
    
    if isempty(ends_5pr)
        ends_5pr(end+1) = max(stops_genes);
    elseif max(ends_5pr)<max(stops_genes)
        % no 5' end before the start of the gene in the buffer region.
        % add an "artificial" 5' end at the beginning of the region of
        % interest.
        ends_5pr(end+1) = max(stops_genes);
    end
    
    if isempty(ends_3pr)
        ends_3pr(end+1) = min(starts_genes);
    elseif min(ends_3pr)>min(starts_genes)
       % if no 3' end after the last gene, generate an artificial one. 
       ends_3pr(end+1) = min(starts_genes);
    end
    
    % get rid of 3' ends before the start of the first gene
    ends_3pr(ends_3pr>max(stops_genes)) = [];
    
    % get rid of 5' ends after the end of the last gene
    ends_5pr(ends_5pr<min(starts_genes)) = [];
    
end


% adding known ends
ends_5pr = sort([ends_5pr known_5pr]);
ends_3pr = sort([ends_3pr known_3pr]);

% remove ends if too close:
min_dist = 10;

diff_ends5 = diff(ends_5pr);
ind_rem5 = find(diff_ends5<min_dist);
ends_5pr(ind_rem5) = [];

diff_ends3= diff(ends_3pr);
ind_rem3 = find(diff_ends3<min_dist);
ends_3pr(ind_rem3) = [];



% removing problematic ends (e.g., processing sites resulting in neighboring 5' and 3' peaks).
bool_rem_5 = true(length(ends_5pr),1);
for i = 1:length(ends_5pr)
    if ismember(ends_5pr(i),x5_rem)
       bool_rem_5(i) = false; 
    end
end
ends_5pr = ends_5pr(bool_rem_5);

bool_rem_3 = true(length(ends_3pr),1);
for i = 1:length(ends_3pr)
    if ismember(ends_3pr(i),x3_rem)
       bool_rem_3(i) = false; 
    end
end
ends_3pr = ends_3pr(bool_rem_3);


% need to flip the data and the positions if on the reverse strand 
if ~strand_region
    x_5 = length(data)-ends_5pr+1;
    x_3 = length(data)-ends_3pr+1;
    data_3pr = fliplr(squeeze(data(1,:,2)));
    data_5pr = fliplr(squeeze(data(1,:,4)));
    x_5 = sort(x_5,'ascend');
    x_3 = sort(x_3,'ascend');

else
   x_5 = ends_5pr;
   x_3 = ends_3pr;
   data_3pr = squeeze(data(1,:,1));
   data_5pr = squeeze(data(1,:,3));
   x_5 = sort(x_5,'ascend');
   x_3 = sort(x_3,'ascend');
end

% feed the data and ends to the isoform quantifier
[levels,termination_probabilities,promoter_strength,...
    isoform_levels, x_sim, simulated_trace] = isoform_quantification_20180130(x_5,x_3,data_3pr,data_5pr,...
    gap_small,gap_read_small,gap_read_large);