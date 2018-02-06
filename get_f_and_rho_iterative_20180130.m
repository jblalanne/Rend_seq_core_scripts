function [f, rho] = get_f_and_rho_iterative_20180130(...
    x,x_5, bool_3, bool_5, levels_subregions, lengths_subregions, max_iter)
% Computes f and rho for region of interest (see definitions below):
% iteratively applies the same method until no more regions are fused
% (rare). 
%
% Inputs:
% x: positions of ends in roi.
% x_5: positions of 5' ends in roi.
% bool_3: bool_3(i) = 1 if x(i) is a 3' end, 0 otherwise.
% bool_5: bool_5(i) = 1 if x(i) is a 5' end, 0 otherwise.
% levels_subregions: mean read counts between ends.
% lengths_subregions: size (nt) of regions between ends.
% max_iter: maximum number of iteration.
%
% Outputs:
% levels_subregions: mean read counts between ends (possibly corrected). 
% f(j): fraction read density past 3' end at x_3(j).
% rho(i): total read density arising from isoforms starting at a 5' end at x_5(i).


% first iteration at 
levels_subregions_0 = levels_subregions;
[levels_subregions_1, f, rho] = get_f_and_rho_20180130(...
    x,x_5,bool_3,bool_5, levels_subregions_0,lengths_subregions);
iter = 2;

% iterate if there are changes (meaning that certain regions were
% effectively fused). 
while (sum(abs(levels_subregions_1-levels_subregions_0)==0)<length(levels_subregions)) ...
        && iter <= max_iter
    [levels_subregions_1, f, rho] = get_f_and_rho_20180130(...
        x,x_5,bool_3,bool_5, levels_subregions_0,lengths_subregions);
    iter = iter+1;
end

if iter==max_iter
   disp('Warning: max number of iterations reached for isoform mapping.'); 
end