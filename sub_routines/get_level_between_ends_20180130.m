function [levels_subregions,lengths_subregions] = ...
    get_level_between_ends_20180130(data_3, x, bool_3, bool_5,...
    gap_small, gap_read_small, gap_read_large)
% gets the mean read counts between reads (using the 3' ends of mapped 
% reads. ensures to leave proper gap next to peaks to avoid possible issues
% with peak shadows.
%
% Inputs:
% data_3: counts per position for 3' ends of mapped reads for roi.
% x: ends in roi.
% bool_3: bool_3(i) = 1 if x(i) is a 3' end, 0 otherwise.
% bool_5: bool_5(i) = 1 if x(i) is a 5' end, 0 otherwise.
% gap_small: gap around peak to be excluded for averaging (~width of peak).
% gap_read_small: smallest read length (serves to set averaging region).
% gap_read_large: smallest read length (serves to set averaging region).
%
% Outputs:
% levels_subregions: mean read counts between ends.
% lengths_subregions: size (nt) of regions between ends.

% 
%% quantify levels between ends
% use the 3' ends of mapped reads for quantification

% initialize
n_ends = length(x);
levels_subregions = zeros(n_ends-1,1);
lengths_subregions = zeros(n_ends-1,1);
bool_empty = false(n_ends-1,1);


% looping through regsion
for i = 1:(n_ends-1)
    
    % the range for averaging depends on the types of ends.
    if bool_3(i)==1 && bool_5(i+1)==1 
        av_range = (x(i)+gap_small):(x(i+1)+gap_read_small);
    elseif bool_3(i)==0 && bool_5(i+1)==1
        av_range = (x(i)+gap_read_large):(x(i+1)+gap_read_small);
    elseif bool_3(i)==1 && bool_5(i+1)==0
        av_range = (x(i)+gap_small):(x(i+1)-gap_small);
    elseif bool_3(i)==0 && bool_5(i+1)==0
        av_range = (x(i)+gap_read_large):(x(i+1)-gap_small);
    end 
    
    % mean level and region sizes
    levels_subregions(i) = mean(data_3(av_range));
    lengths_subregions(i) = length(av_range);
    if lengths_subregions(i)==0
           bool_empty(i)=1; 
    end
end


% for 3' ends too close to 5 ends: replace by mean of neighboring downstream level
ind_empty = find(bool_empty);
for i = 1:length(ind_empty)
    if ind_empty(i) ==1
        levels_subregions(ind_empty(i)) =levels_subregions(ind_empty(i)+1);
    elseif ind_empty(i) == n_ends-1
        levels_subregions(ind_empty(i)) = levels_subregions(ind_empty(i)-1);
    else
        levels_subregions(ind_empty(i)) = 0.5*(levels_subregions(ind_empty(i)-1)+levels_subregions(ind_empty(i)+1));
    end
end

% 
% 
% levels = zeros(n_ends-1,1);
% level_sizes = zeros(n_ends-1,1);
% bool_empty = false(n_ends-1,1);
% for i = 1:(n_ends-1)
%     
%     if bool_3(i)==1 && bool_5(i+1)==1
%         levels(i) = mean(data_3((x(i)+gap_small):(x(i+1)+gap_read_small)));
%         level_sizes(i) = length((x(i)+gap_small):(x(i+1)+gap_read_small));
%     elseif bool_3(i)==0 && bool_5(i+1)==1
% %         levels(i) = mean(data_3((x(i)+gap_read_large):(x(i+1)+gap_read_small)));
% try
%         levels(i) = mean(data_5((x(i)+gap_small):(x(i+1)-gap_small)));
%         level_sizes(i) = length((x(i)+gap_read_large):(x(i+1)+gap_read_small));
% catch
%     bla
%     
% end
%     elseif bool_3(i)==1 && bool_5(i+1)==0
%         levels(i) = mean(data_3((x(i)+gap_small):(x(i+1)-gap_small)));
%         level_sizes(i) = length((x(i)+gap_small):(x(i+1)-gap_small));
%     elseif bool_3(i)==0 && bool_5(i+1)==0
%         levels(i) = mean(data_3((x(i)+gap_read_large):(x(i+1)-gap_small)));
%         level_sizes(i) = length((x(i)+gap_read_large):(x(i+1)-gap_small));
%         if level_sizes(i)==0
%            bool_empty(i)=1; 
%         end
%     end 
% end
% 
% 
% 
% % dealing with terminators too close to promoter: replace by nearest
% % downstream level
% ind_empty = find(bool_empty);
% for i = 1:length(ind_empty)
%    levels(ind_empty(i)) = levels(ind_empty(i)+1); 
% end
% 
