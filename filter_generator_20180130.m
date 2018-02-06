function filter = filter_generator_20180130(length_data,half_width,filter_type,gap)
% generates filter for convolution (for averaging).
%
% Inputs:
% length_data: length of data matrix to be average (e.g., chromosome size).
% half_width: size of filter for averaging.
% filter_type: determines the type of filter.
% gap: gap variable for pitted filters.
%
% Output:
% filter: filter for averaging by FFT convolution. 


% initialize filter.
filter = zeros(length_data,1);

% generates filter 
if strcmp(filter_type,'central even')
    filter(1+(ceil(length_data/2)-half_width:ceil(length_data/2)+half_width)) = ...
        1/(2*half_width+1)*ones(2*half_width+1,1);
    
elseif strcmp(filter_type,'central')
    filter(1+(ceil(length_data/2)-half_width:ceil(length_data/2)+half_width)) = ...
        1/(2*half_width+1)*ones(2*half_width+1,1);
    
elseif strcmp(filter_type,'pitted central')
    filter(1+(length_data/2-half_width-gap:length_data/2-1-gap)) = ...
        1/(2*half_width)*ones(half_width,1);
    
    filter(1+(length_data/2+1+gap:length_data/2+half_width+gap)) = ...
        1/(2*half_width)*ones(half_width,1);
    
elseif strcmp(filter_type,'forward')
    filter(length_data/2+1:length_data/2+2*half_width) = 1/(2*half_width);
    
elseif strcmp(filter_type,'reverse')
    filter(length_data/2-2*half_width:length_data/2-1) = 1/(2*half_width);
else
    disp('Filter type not recognized.');
    
end

end