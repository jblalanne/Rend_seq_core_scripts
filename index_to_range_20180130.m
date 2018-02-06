function range = index_to_range_20180130(x,dx_long,dx_short,index,genome_size)
% determines the averaging range given the type of end (upstream of 3' ends,
% downstream of 5' ends) and strand.
%
% Input: 
% x: position.
% dx_long: distance to end of range.
% dx_short: distance to start of range.
% index: determines the type of end (1: 3'f, 2: 3'r, 3: 5'f, 4:5'r).
% genome_size: size of genome (to avoid going beyond). 
%
% Output:
% range: range of interest.

if index==1             % 3f
    range = (x-dx_long):(x-dx_short);
elseif index == 2       % 3r
    range = (x+dx_short):(x+dx_long);
elseif index == 3       % 5f
    range = (x+dx_short):(x+dx_long);
elseif index == 4       % 5r
    range = (x-dx_long):(x-dx_short);
end

% dealing with positions out of bounds (negative or larger
% than genome size).
if min(range)<0
    range = 1:range(end);
elseif max(range)>genome_size
    range = range(1):genome_size;
end