function isoform_levels = get_isoform_levels_20180130(x,x_5,x_3,f,rho)
% compute the level of mRNA isoforms from f and rho (see Methods of
% publication). 
%
% Inputs:
% x: positions of ends in roi.
% x_5: positions of 5' ends in roi.
% x_3: positions of 3' ends in roi.
% f(j): fraction read density past 3' end at x_3(j).
% rho(i): total read density arising from isoforms starting at a 5' end at x_5(i).
%
% Output
% isoform_levels: isoform level array, entry at (i,j) equals abundance of
%                 mRNA isoform starting at x_5(i) and ending at x_3(j).


% initialization
isoform_levels = zeros(length(x_5),length(x_3));
    
% looping through 5' ends
for i = 1:length(x_5)
    
    % index of end of interest
    ind5 = find(x==x_5(i));
    
    % indices of 5' ends upstream of x_5(i).
    ind5_up = 1:(i-1);
    
    % indices of 3' ends downstream of x_5(i)
    ind3_down = find(x_3>x_5(i));
   
    % applying 
    for j = 1:length(ind3_down)
        non_term_factor = 1;
        for k = 1:(j-1)
            non_term_factor = non_term_factor*(1-f(ind3_down(k)));
        end
        isoform_levels(i,ind3_down(j)) = rho(i)*non_term_factor*f(ind3_down(j));
    end
end
