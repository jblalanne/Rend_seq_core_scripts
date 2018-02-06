function [data_3,data_5, x_3, x_5] =...
    get_ends_roi_20180130(start_region,stop_region,strand_region,...
    z_thresh,width_z,z_peak_oi,data_oi,known_5,known_3,spurious_5,spurious_3)
% Gets the ends (identified from peak z score) in region of interest.
% 
% Input:
%


% defining region of interest
range = start_region:stop_region;

% looking for 3' and 5' peaks
if strand_region
    ends_3 = find(z_peak_oi(:,1)>z_thresh);
    ends_5 = find(z_peak_oi(:,3)>z_thresh);
else
    ends_3 = find(z_peak_oi(:,2)>z_thresh);
    ends_5 = find(z_peak_oi(:,4)>z_thresh);    
end

% keep only the peaks in roi
ends_3 = ends_3(ends_3>min(range) & ends_3<max(range));
ends_5 = ends_5(ends_5>min(range) & ends_5<max(range));


% adding known ends
ends_5 = sort([ends_5 known_5]);
ends_3 = sort([ends_3 known_3]);

% remove redundant positions
temp = ends_3;
temp(diff(temp)<=width_z) = [];
ends_3 = temp;
temp = ends_5;
temp(diff(temp)<=width_z) = [];
ends_5 = temp;


% removing spurious ends (e.g., processing sites resulting in immediate 
% neighboring 5' and 3' peaks missed by automatic peak finding). 
bool_rem_5 = true(length(ends_5),1);
for i = 1:length(ends_5)
    if ismember(ends_5(i),spurious_5)
       bool_rem_5(i) = false; 
    end
end
ends_5 = ends_5(bool_rem_5);

bool_rem_3 = true(length(ends_3),1);
for i = 1:length(ends_3)
    if ismember(ends_3(i),spurious_3)
       bool_rem_3(i) = false; 
    end
end
ends_3 = ends_3(bool_rem_3);


% if first end is a 3' end: artificially set x=1 as a 5' end.
if strand_region
    if min(ends_3)<min(ends_5)
        ends_5 = [min(range) ends_5];
    end
else
    if max(ends_3)>max(ends_5)
        ends_5 = [max(range) ends_5];
    end
end

% if last end is a 5' end, artificially set last range position as a 3' end.
if strand_region
    if max(ends_3)<max(ends_5)
        ends_3 = [ends_3; max(range)]; 
    end
else
    if min(ends_3)>min(ends_5)
        ends_3 = [min(range); ends_5];
    end
end


% flip the data and the positions if on the reverse strand 
if ~strand_region
    x_5 = length(data_oi)-ends_5+1;
    x_3 = length(data_oi)-ends_3+1;
    data_3 = flipud(squeeze(data_oi(:,2)));
    data_5 = flipud(squeeze(data_oi(:,4)));
    x_5 = sort(x_5,'ascend');
    x_3 = sort(x_3,'ascend');

else
   x_5 = ends_5;
   x_3 = ends_3;
   data_3 = squeeze(data_oi(:,1));
   data_5 = squeeze(data_oi(:,3));
   x_5 = sort(x_5,'ascend');
   x_3 = sort(x_3,'ascend');
end

