function [average_data, z_peak] = compute_peak_z_score_20180130(data,half_width,...
    average_threshold,filter_type,gap)
% computes the modified peak z score of a 1-D signal.
%
% Inputs:
% data: data variable containing read counts per position.
% half_width: half-width of filter for averaging.
% average_threshold: threshold in average read count for consideration.
% filter_type: type of filter for z score computation.
% 
% Outputs:
% average_data: average data (with same filter as z score computation).
% z_score: modified peak z score for data. 


% obtain filter for convolution.
filter = filter_generator_20180130(length(data),half_width,filter_type,gap);

% discrete Fourier transform of the data and of the filter.
ft_data = fft(data);
ft_filter = fft(filter);

% discrete Fourier transform of the square of the data.
ft_data_squared = fft(data.^2);

% perform convolutions via fast Fourier transform.
average_data = fftshift(ifft(ft_filter'.*ft_data));
std_data = sqrt(fftshift(ifft(ft_filter'.*ft_data_squared)) - ...
    average_data.^2);

% ensure that if the average is 0, you avoid NAN
z_peak = (average_data>average_threshold).*(data-average_data)./(std_data + (average_data==0));


