
% The current script illustrates a simple analysis work flow for Rend-seq
% data using the B. subtilis LB, 25 s fragmentation, dataset. The dataset
% used below is included in the GitHub submission.


%% import .wig and compute z score
% Import data and compute peak and step z score.

% parameters for peak z score and step z score calculation
half_width_z_score = 50;    % half width of window for averaging for peak z score
half_width_step = 100;      % half width of window for averaging for step z score
average_threshold = 0.25;   % average read density (read/nt) threshold for consideration
gap_z = 2;                  % gap left out (both sides) of central position for peak z score 
gap_step = 3;               % gap left out (both sides) of central position for step z score 

% full path of the directory containing the wig files. Needs to be modified
% to directory of interest
data_dir = 'my_dir/example_data';

% file name for different data series: wild card takes care of 3f, 3r, 5f, 5r (assumes that order).  
% (more can be added: files{2} = 'mydata_*.wig')
files{1} = 'Bacillus_subtilis_WT_Rend_seq_LB_25s_frag_pooled_*.wig';

% genome size in nucleotides (required for array initialization, here for B. subtilis):
genome_size = 4215606;

% get_data_compute_statistic_simple: imports the wig file in a Matlab array
% and computes the peak and z score. Takes about ~30 s per dataset.
[data,z_peak,z_step] =  get_data_compute_statistic_20180130(...
    files,data_dir,...
    half_width_z_score,half_width_step,...
    average_threshold,gap_step,gap_z,genome_size);


%% peak calling

% parameters for peak calling
z_thresh = 12;      % positions with peak z score>z_thresh considered "peaks"
d_thresh = 3;       % positions within d_thresh aggregated to highest neighboring z score position.

% getting the peaks positions
peak_positions = get_z_peaks_20180130(z_peak,z_thresh,d_thresh);


%% generating the end-enrichment file for shadow removal

l_average = 50;     % size of region for averaving for end-enrichment calculation.
width_around_peak = 10;
ee_name_addition = 'end_enrichment_shadow_removal.txt';    % suffix to file name of end-enrichment file for shadow removal.

% prints end-enrichment file: 
%       column 1: type of ends, column 2: position, column 3: end-enrichment
print_end_enrichment_shadow_removal_20180130(files,ee_name_addition,...
    peak_positions,data,width_around_peak,l_average,genome_size);

% the wigs without shadow can then be generated using the python script by:
% python align_to_wig_no_shadow_20180130.py bowtie_align_file new_wig_name end_enrichment_file species 

