function counts = extract_wig_20180130(files_name,data_dir,genome_size)
% reads wig files.
%
% Input:
% files_name: name of Rend-seq wig files (wild card for 3f, 3r, 5f, 5r).
% data_dir: full path of directory containing wig files.
% genome_size: size of genome/chromosome of interest.
%
% Output:
% counts: 3D array with 
        % dimension 1: different data series (same size as files).
        % dimension 2: read counts per position (genome size).
        % dimension 3: type of ends (in order):
        %                   3' end of mapped reads forward,
        %                   3' end of mapped reads reverse,
        %                   5' end of mapped reads forward,
        %                   5' end of mapped reads reverse.


% current directory
cwd = pwd;

% moves to directory of interest
cd(data_dir);
files = dir(files_name);

if isempty(files)
    disp('Files not found.');
    counts = NaN;
else

    % initializing array
    counts = zeros(length(files),genome_size);
    
    % reading data & converting wig to data array where index is genome
    % position
    for i = 1:length(files)
        fid = fopen(files(i).name);
        disp(files(i).name)
        
        content = textscan(fid,'%d %f','Delimiter','\t','Headerlines',2);
        indices = content{1};
        reads = content{2};
        counts(i,indices) = reads;
        
        fclose(fid);
    end
end

cd(cwd)