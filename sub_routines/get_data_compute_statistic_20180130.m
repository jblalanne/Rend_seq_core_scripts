function [data,z_peak,z_step,average_data] = ...
    get_data_compute_statistic_20180130(files,data_dir,half_width_z_score,half_width_step,...
    average_threshold,gap_step,gap_z,genome_size)
% Extracts wig files and computes peak and step z scores. 
%
% % % Inputs:
% files: cell array with wig file names (one per dataset with wild-card).
% data_dir: full path of directory containing wig files.
% half_width_z_score: half width of window for averaging for peak z score
% half_width_step = 100: half width of window for averaging for step z score
% average_threshold: average read density (read/nt) threshold for consideration
% gap_z: gap left out (both sides) of central position for peak z score 
% gap_step: gap left out (both sides) of central position for step z score 
%
% Outputs:
% data: 3D data array containing the read counts per position: 
    % dimension 1: data series
    % dimension 2: genome coordinate
    % dimension 3: type of signal (in order): 3' end mapped reads forward,
    %              3' end mapped reads reverse,
    %              5' end mapped reads forward,
    %              5' end mapped reads reverse.
% z_stat: same as data (same organization) with peak z score.
% step_stat: same as data (same organization) with step z score.
% average: same as data (same organization) with average read count.


% % % % % % % % % % % % % % %
% % % Importing the wig % % %
% % % % % % % % % % % % % % %

% initialization
data = zeros(length(files),genome_size,4);
% index key for 3rd dimension:
% 1: 3' forward
% 2: 3' reverse
% 3: 5' forward
% 4: 5' reverse

for i = 1:length(files)
    counts = extract_wig_20180130(files{i},data_dir,genome_size);
    data(i,:,1) = counts(1,:);
    data(i,:,2) = counts(2,:);
    data(i,:,3) = counts(3,:);
    data(i,:,4) = counts(4,:);
end

% requires even array for FFT
if mod(length(data),2)==1
    data(:,end+1,:) = 0;
end


% % % % % % % % % % % % % % % % %
% % % Computing statistics  % % %
% % % % % % % % % % % % % % % % %

% initialization
z_peak = zeros(length(files),size(data,2),4);
average_data = zeros(length(files),size(data,2),4);
z_step = zeros(length(files),size(data,2),4);

for i = 1:length(files)
    for j = 1:4
        % computation of peak and step z score.
        [average_data(i,:,j),z_peak(i,:,j)] = compute_peak_z_score_20180130(squeeze(data(i,:,j)),half_width_z_score,...
            average_threshold,'pitted central',gap_z);
        z_step(i,:,j) = compute_step_z_score_20180130(squeeze(data(i,:,j)),...
            half_width_step,average_threshold,gap_step);
    end
end