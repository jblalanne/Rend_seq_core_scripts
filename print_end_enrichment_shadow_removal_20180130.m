function print_end_enrichment_shadow_removal_20180130(files,ee_name_addition,...
    peak_positions,data,width_around_peak,l_average,genome_size)
% Takes peak positions and computes end-enrichment at each position. Prints
% file used for shadow removal. Not aware of additional neighboring peaks. 
%
% Input: 
% files_name: name of dataset
% ee_name_addition: name addition for end-enrichment file. 
% peak_positions: data structure containing the position of peaks for each
%                 dataset
% width_around_peak: region surrounding peaks for peak value determination.
% l_average: range for averaging. 
% genome_size: size of genome. 
%
% Output: 
% no variable output.
% Prints to current directory the end-enrichment file for shadow removal,
% which contains three columns of information. column 1: type of ends, 
% column 2: position, column 3: end-enrichment.

for i = 1:length(files)
    
    % name of end-enrichment file
    file_name_end_enrichment = [strrep(files{i},'*.wig','') ee_name_addition];
    
    % opening file
    fid = fopen(file_name_end_enrichment,'w');
    
    for j = 1:4
        for k = 1:length(peak_positions{i,j})

            % determining peak width from sum of reads above half-max.
            peak_region = (peak_positions{i,j}(k)-width_around_peak):(peak_positions{i,j}(k)+width_around_peak);
            if min(peak_region)<0
                peak_region = 1:peak_region(end);
            elseif max(peak_region)>genome_size
                peak_region = peak_region(1):genome_size;
            end
            
            data_peak = squeeze(data(i,peak_region,j));
            peak_value = sum(data_peak(data_peak>0.5*max(data_peak)));
            
            % range: upstream for 3' ends, downstream for 5' ends
            range_surrounding_peak = index_to_range_20180130(peak_positions{i,j}(k),l_average+2*width_around_peak,2*width_around_peak,j,genome_size);
            mean_surrounding_peak = mean(data(i,range_surrounding_peak,j));
            enrichment_factor = peak_value/mean_surrounding_peak;
            
            % printing the end-enrichment for each peak. 
            fprintf(fid,'%d\t%8d\t%.2f\n',j,peak_positions{i,j}(k),enrichment_factor); 
        end
    end
    fclose(fid);
    fprintf(sprintf('Done obtaining end-enrichment for %s.\n',strrep(files{i},'*.wig','')));
    
end


