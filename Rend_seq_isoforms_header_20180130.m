
% The current script illustrates a simple analysis workflow for Rend-seq
% data using the B. subtilis LB, 25 s fragmentation, dataset. The dataset
% below is included in the GitHub submission.


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
files{1} = 'Bacillus_subtilis_WT_Rend_seq_LB_25s_frag_no_shadow_pooled_*.wig';

% genome size in nucleotides (required for array initialization, here for B. subtilis):
genome_size = 4215606;

% get_data_compute_statistic_simple: imports the wig file in a Matlab array
% and computes the peak and z score. Takes about ~30 s per dataset.
[data,z_peak,z_step] =  get_data_compute_statistic_20180130(...
    files,data_dir,...
    half_width_z_score,half_width_step,...
    average_threshold,gap_step,gap_z,genome_size);



%% Example display of Rend-seq data with mRNA isoform quantification

% Enter the start and end coordinates of region of interest, as well as the
% strand. Run cell. 

% example: rpsP operon.
start_region = 1671780;
stop_region = 1676420;
strand_region = 1;          

% % example 2 (uncomment and run for display): qoxA operon (note that the definition is such that start<stop,
% % even on the reverse strand). 
% start_region = 3914268;
% stop_region = 3918624;
% strand_region = 0;

% in the case of multiple data series, choose which one to display with
% ind_oi.
ind_oi = 1;
z_peak_oi = squeeze(z_peak(ind_oi,:,:));
data_oi = squeeze(data(ind_oi,:,:));
z_thresh = 12;
width_z = 3;

% Genbank annotation file and directory for gene annotation. Directory
% needs to be modified. 
annotation_file = 'NC_000964.3.faa';
annotation_dir = 'my_dir/example_data';

% if there are missed peaks: enter manually (from position on chromosome)
known_5 = [];
known_3 = [];

% if there are spurious peaks to be removed: enter manually (from position on chromosome)
spurious_5 = [];
spurious_3 = [];


% identification, quantification and plot of region of interest. 
[x_3, x_5, isoform_levels] = quantify_plot_isoforms_20180130(start_region,stop_region,...
    strand_region,z_thresh,width_z,z_peak_oi,data_oi,known_5,known_3,...
    spurious_5,spurious_3,genome_size,annotation_file,annotation_dir);
% outputs of above:
% x_3: position of 3' ends in region of interest.
% x_5: position of 5' ends in region of interest.
% isoform_levels: isoform level array, entry at (i,j) equals abundance of
%                 mRNA isoform starting at x_5(i) and ending at x_3(j).

