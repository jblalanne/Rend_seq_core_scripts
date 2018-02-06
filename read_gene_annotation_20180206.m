function [start_f, stop_f, genes_f, start_r, stop_r, genes_r] =...
    read_gene_annotation_20180206(annotation_file,dir_name)
% Reads Genbank annotation file and returns gene names, strand, start and
% stop. Start position is always smaller than stop position, even on the
% reverse strand (meaning that the start position corresponds to the stop
% codon on the reverse strand, and vice versa).
%
% Inputs:
% file_name: name of Genbank .faa file.
% dir_name: full path of directory containing .faa file.
%
% Outputs:
% start_f: start positions (in nt on the chromosome) of genes on
%               forward strand. 
% stop_f: stop positions (in nt on the chromosome) of genes on
%               forward strand. 
% genes_f: name of genes on forward strand.
% start_r: start positions (in nt on the chromosome) of genes on
%               reverse strand. Correspond to ends of gene.
% stop_r: stop positions (in nt on the chromosome) of genes on
%               forward strand. Correspond to starts of gene.
% genes_r: name of genes on reverse strand.


% line by line reading 
gene_counter = 0;
strand = [];
gene_name = [];
start = [];
stop = [];

curr_dir = pwd;
cd(dir_name);
fid = fopen(annotation_file);
tline = fgetl(fid);


while ischar(tline) && ~isempty(tline)
%     disp(tline)

    if strcmp(tline(1),'>')
        gene_counter = gene_counter+1;
        
        % getting the gene name
        startIndex = regexp(tline,'gene=')+5;
        stopIndex = regexp(tline,'protein=')-4;
        
        gene_info = tline(startIndex:stopIndex);
        
        % if multiple gene names, pick the first one.
        final_stopIndex = regexp(gene_info,',')-1;
        
        if isempty(final_stopIndex)
            gene_name{gene_counter} = gene_info;
        elseif isempty(gene_info)
            gene_name{gene_counter} = gene_info(1:final_stopIndex);
        else
            gene_name{gene_counter} = sprintf('gene%d',gene_counter);
            i
            break
        end
      
        try
        % getting the strand and position
        startIndex = regexp(tline,'location=complement(')+20;
        if ~isempty(startIndex)
            strand(gene_counter) = 0;
            position_string = tline(startIndex:end);
            position_string(position_string=='.')='-';
           
            ind1 = regexp(position_string,'-->')-1;
            ind1b = regexp(position_string,'--')-1;
            ind2 = regexp(position_string,')')-1;
            
            

            if isempty(ind1)
                start(gene_counter) = str2num(position_string(1:ind1b));
                stop(gene_counter) = str2num(position_string(ind1b+3:ind2));
            else
                start(gene_counter) = str2num(position_string(1:ind1));
                stop(gene_counter) = str2num(position_string(ind1+4:ind2));
            end
            
        else
            strand(gene_counter) = 1;
            startIndex = regexp(tline,'location=')+9;
            position_string = tline(startIndex:end);
            position_string(position_string=='.')='-';

            ind1 = regexp(position_string,'-->')-1;
            ind1b = regexp(position_string,'--')-1;
            ind2 = regexp(position_string,']')-1;

            if isempty(ind1)
                start(gene_counter) = str2num(position_string(1:ind1b));
                stop(gene_counter) = str2num(position_string(ind1b+3:ind2));
            else
                start(gene_counter) = str2num(position_string(1:ind1));
                stop(gene_counter) = str2num(position_string(ind1+4:ind2));
            end
            
            
        end
        catch
            disp('Problem:');
            disp(tline);
        end
    end
    tline = fgetl(fid);
    
end

fclose(fid);
cd(curr_dir)


genes_f = gene_name(strand==1)';
genes_r = gene_name(strand==0)';

start_f = start(strand==1);
start_r = start(strand==0);

stop_f = stop(strand==1);
stop_r = stop(strand==0);

