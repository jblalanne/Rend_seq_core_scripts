function [x,bool_3,bool_5] = reorder_positions_20180130(x_3,x_5)
% orders 3' and 5' ends based on position, irrespectivd of type
%
% Input:
% x_3: positions of 3' ends.
% x_5: positions of 5' ends.
%
% Output:
% x: positions of all ends.
% bool_3: bool_3(i) = 1 if x(i) is a 3' end, 0 otherwise.
% bool_5: bool_5(i) = 1 if x(i) is a 5' end, 0 otherwise.



x = [x_3' x_5'];
[x,ind_sort] = sort(x,'ascend');

bool_3 = [ones(size(x_3')) zeros(size(x_5'))];
bool_3 = bool_3(ind_sort);
bool_5 = ~bool_3;




% % temporary array progressively depleted
% temp_x3 = x_3;
% temp_x5 = x_5;
% 
% % number of ends
% n_ends = (length(x_3)+length(x_5));
% 
% % initializting
% x = zeros(n_ends,1);
% bool_5 = zeros(n_ends,1);
% bool_3 = zeros(n_ends,1);
% 
% 
% % looping through and ordering. 
% for i = 1:n_ends
%     
%     if ~isempty(temp_x5) && isempty(temp_x3)
%         if temp_x3(1)<temp_x5(1)
%             x(i) = temp_x3(1);
%             bool_3(i) = 1;
%             temp_x3 = temp_x3(2:end);
%         elseif temp_x3(1)>temp_x5(1)
%             x(i) = temp_x5(1);
%             bool_5(i) = 1;
%             temp_x5 = temp_x5(2:end);
%         end
%     else
%         
%         x(i)=temp_x3(1);
%         bool_3(i) = 1;
%         temp_x3 = temp_x3(2:end);
%     end
% end
