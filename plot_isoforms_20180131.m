function plot_isoforms_20180131(x_rec,rec_trace,x_5,x_3,data_5,data_3,...
    strand_region,isoform_levels,genome_size,annotation_file,annotation_dir)
% Displays Rend-seq data, with gene positions and identified mRNA isoforms
% with corresponding abundance (as determined by Rend-seq) for a specified
% region of interest.
%
% Inputs:
% x_rec: positions for reconstructed Rend-seq trace, see
%       isoform_quantification_20180130.m.
% rec_trace: reconstructed Rend-seq trace, see isoform_quantification_20180130.m.
% x_5: identified 5' ends in region of interest.
% x_3: identified 3' ends in region of interest.
% data_5: 5' ends of mapped reads for strand of interest.
% data_3: 3' ends of mapped reads for strand of interest.
% strand_region: strand of interest.
% isoform_levels: isoform abundance, see isoform_quantification_20180130.m.
% genome_size: size of genome (in nt).
% annotation_file: file name for GenBank annotation file.
% annotation_dir: full path of directory containing annotation file. 


%% display Rend-seq data
Y_max = 2.5E4;      % maximum of read count range

figure;
subplot(9,1,1:6);
hold on

extra_dist = 100;

y_lines = logspace(0,4,5);
for i = 1:length(y_lines)
   plot([x_rec(1)-extra_dist x_rec(end)+extra_dist],y_lines(i)*[1 1],'-','Color',[0.8 0.8 0.8]) 
   text(0,log10(y_lines(i))*0.23,sprintf('10^%d',log10(y_lines(i))),...
       'HorizontalAlignment','right','FontSize',12,'Color',[0 0 0],'Units','Normalized')
end

text(-0.06,0.5,'Read counts',...
    'HorizontalAlignment','center','FontSize',15,'Color',[0 0 0],'Rotation',90,'Units','Normalized')

x_range = (min(x_rec)-extra_dist):(max(x_rec)+extra_dist);

color_3 = [0 122 255]/255;
color_5 = [255 137 0]/255;
stairs(x_range,data_5(x_range),'-','Color',color_5,'LineWidth',1.5)
stairs(x_range,data_3(x_range),'-','Color',color_3,'LineWidth',1.5)
stairs(x_rec,rec_trace,'-k','LineWidth',1.5)
set(gca,'XLim',[min(x_range) max(x_range)],'YLim',[1 Y_max])
set(gca,'YScale','log')
axis off


%% display gene positions in region of interest

% obtaining gene positions from annotation file.
[start_f, stop_f, genes_f, start_r, stop_r, genes_r] = read_gene_annotation_20180206(annotation_file,annotation_dir);


% find genes in the roi. 
if strand_region
    ind_genes = find( (start_f>min(x_rec)) & ...
        (stop_f<max(x_rec)));
else
    min_gene_other = genome_size-max(x_rec);
    max_gene_other = genome_size-min(x_rec);

    ind_genes = find(start_r>min_gene_other & ...
        stop_r<max_gene_other );
end


% changing coordinates if on the reverse strand (flipping from left to
% right). 
if ~strand_region
    starts_roi = genome_size-stop_r(ind_genes)+1;
    stops_roi = genome_size-start_r(ind_genes)+1;
    genes_roi = {genes_r{ind_genes}};
else
    starts_roi = start_f(ind_genes);
    stops_roi = stop_f(ind_genes);
    genes_roi = {genes_f{ind_genes}};
end


% gene arrow display parameters
w = 10;
dw = 2;
l = 100;

subplot(9,1,7);
hold on
for i = 1:length(ind_genes)
    [x,y] = arrow_function_20180131(starts_roi(i),stops_roi(i),w,dw,l);
    patch(x,y,'w')
    text(mean([starts_roi(i) stops_roi(i)]),w/2,genes_roi{i},...
        'HorizontalAlignment','center','Color','k','FontSize',12,'FontAngle','italic')
end

set(gca,'XLim',[min(x_range) max(x_range)],'YLim',[-0.5*w 1.1*w+dw])
axis off


%% display identified isoform with quantification

% color scale for isoforms
isoform_levels_norm = isoform_levels/max(max(isoform_levels));
isoform_levels_norm(isoform_levels_norm<0) = 0;

% isoform with minimal relative abundance displayed is 10^-3.
min_color = -3;
min_level = log10(min(min(isoform_levels_norm(isoform_levels_norm>0))));
min_level_above_cut = log10(min(min(isoform_levels_norm(isoform_levels_norm>10^min_color))));
min_color = max([min_level min_level_above_cut]);
max_color = 0;


% color scale for isoform abundance
n_colors = 500;
color_index = zeros(size(isoform_levels));
for i = 1:length(x_5)
    for j = 1:length(x_3)
        if log10(isoform_levels_norm(i,j))<=min_color
            color_index(i,j) = 1;
        elseif log10(isoform_levels_norm(i,j))>=max_color
            color_index(i,j)=n_colors;
        else
            color_index(i,j) = ceil(n_colors*(log10(isoform_levels_norm(i,j))-min_color)/(max_color-min_color))+1;
        end
    end
end
c = [ ones(n_colors,1) zeros(n_colors,2) ];
for i = 1:n_colors
    c(i,:) = pale_RGB_20180131(c(i,:),1-i/n_colors);
end


% plotting isoforms with abundance
subplot(9,1,8:9);
hold on
counter = 1;
width = 6.5;
for i = 1:length(x_5)
    for j = 1:length(x_3)
        if log10(isoform_levels_norm(i,j))>=min_color
            [x,y] = rectangle_points_20180131(x_5(i),x_3(j),-counter*10,width);
            patch(x,y,c(color_index(i,j),:),'LineWidth',0.5)
            text(x_3(j)+25,-counter*10-2,sprintf('%.3f',isoform_levels_norm(i,j)),'FontSize',8,...
                'VerticalAlignment','baseline')
            counter=counter+1;
        end
    end
end
colormap(c)
caxis([min_color max_color])
set(gca,'XLim',[min(x_range) max(x_range)])
set(gca,'YLim',[-counter*10 0])
axis off
