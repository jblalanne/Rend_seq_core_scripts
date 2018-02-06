function [levels_subregions, f, rho] = get_f_and_rho_20180130(...
    x,x_5,bool_3,bool_5, levels_subregions,lengths_subregions)
% Computes f and rho for region of interest (see definitions below).
%
% Inputs:
% x: positions of ends in roi.
% x_5: positions of 5' ends in roi.
% bool_3: bool_3(i) = 1 if x(i) is a 3' end, 0 otherwise.
% bool_5: bool_5(i) = 1 if x(i) is a 5' end, 0 otherwise.
% levels_subregions: mean read counts between ends.
% lengths_subregions: size (nt) of regions between ends.
%
% Outputs:
% levels_subregions: mean read counts between ends (possibly corrected). 
% f(j): fraction read density past 3' end at x_3(j).
% rho(i): total read density arising from isoforms starting at a 5' end at x_5(i).



% number of ends in roi
n_ends = length(x);

%% f: fraction of read density past identified 3' ends
f = [];
for i = 1:n_ends-1
    if bool_3(i)
        f(end+1) = 1-(levels_subregions(i)/levels_subregions(i-1));
        
        % there is negative f: fuse neighboring regions.
        if f(end)<0
            f(end) = 0;
            average_level = (levels_subregions(i)*lengths_subregions(i) + levels_subregions(i-1)*lengths_subregions(i-1))/(lengths_subregions(i)+lengths_subregions(i-1));
            levels_subregions(i) = average_level;
            levels_subregions(i-1) = average_level;
        end
    end
end
% assume last end terminates 100%.
f(end+1) = 1;


%% rho: get total read density increase past 5' end

% initialization
rho = zeros(length(x_5),1);

% looping through x_5
for i = 1:length(x_5)
    % index of end of interest
    ind5 = find(x==x_5(i));
    if ind5==1
        rho(i) = levels_subregions(ind5); % if x_5(i) is the first 5' end
    else
        rho(i) = levels_subregions(ind5)-levels_subregions(ind5-1);
    end

end

% negative rho: fuse neighboring regions.
ind_neg_prom = find(rho<0);
ind_5 = find(bool_5);
for j = 1:length(ind_neg_prom)
    i = ind_5(ind_neg_prom(j));
    average_level = (levels_subregions(i)*lengths_subregions(i) + levels_subregions(i-1)*lengths_subregions(i-1))/(lengths_subregions(i)+lengths_subregions(i-1));
    levels_subregions(i) = average_level;
    levels_subregions(i-1) = average_level;
end

