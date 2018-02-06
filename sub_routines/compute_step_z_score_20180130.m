function z_step = compute_step_z_score_20180130(data,half_width,av_thresh,gap)
% computes the step z score of a 1-D signal.
%
% Inputs:
% data: data variable containing read counts per position.
% half_width: half-width of filter for averaging.
% average_threshold: threshold in average read count for consideration.
% 
% Outputs:
% average_data: average data (with same filter as z score computation).
% z_score: modified peak z score for data. 


length_data = length(data);

% filters
right_filter = zeros(length_data,1);
right_filter(length_data/2+1+gap:length_data/2+half_width+gap) = ...
    1/(half_width)*ones(half_width,1);

left_filter = zeros(length_data,1);
left_filter((length_data/2-half_width-gap:length_data/2-1-gap)) = ...
    1/(half_width)*ones(half_width,1);

% computing averages by convolution
right_av = convolve_20180130(data,right_filter);
right_av_2 = convolve_20180130(data.^2,right_filter);
right_var = right_av_2-right_av.^2;

left_av = convolve_20180130(data,left_filter);
left_av_2 = convolve_20180130(data.^2,left_filter);
left_var = left_av_2-left_av.^2;

% want the counts to be above some threshold on at least one side
average_bool = right_av>av_thresh | left_av>av_thresh;

z_step = average_bool.*(right_av-left_av)./sqrt(left_var + right_var+~average_bool);

end