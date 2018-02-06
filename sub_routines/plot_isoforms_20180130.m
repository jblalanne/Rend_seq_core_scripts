function plot_isoforms_20180130(x_rec,rec_trace,data_5,data_3,strand_region)



%% data with 
Y_max = 2.5E4;
buffer = 25;

subplot(9,1,1:6);
hold on

extra_dist = 100;

y_lines = logspace(0,4,5);
for i = 1:length(y_lines)
   plot([x_rec(1)-extra_dist x_rec(end)+extra_dist],y_lines(i)*[1 1],'-','Color',[0.8 0.8 0.8]) 
   text(0,log10(y_lines(i))*0.2025+0.15,sprintf('10^%d',log10(y_lines(i))),...
       'HorizontalAlignment','right','FontSize',12,'Color',[0 0 0],'Units','Normalized')
end

text(-0.06,0.5,'Read counts',...
    'HorizontalAlignment','center','FontSize',15,'Color',[0 0 0],'Rotation',90,'Units','Normalized')


color_3 = [0 114 178]/255;
color_5 = [213 94 0]/255;
stairs(x_rec,data_5(x_rec),'-','Color',pale_RGB(color_5,0.5),'LineWidth',1.5)
stairs(x_rec,data_3(x_rec),'-','Color',pale_RGB(color_3,0.5),'LineWidth',1.5)
stairs(x_rec,rec_trace,'-k','LineWidth',1.5)
set(gca,'XLim',[(min(x_rec)-buffer) (max(x_rec)+buffer)],'YLim',[0.25 Y_max])
set(gca,'YScale','log')
axis off


%% genes in operon

[start_forward, stop_forward, genes_f, ...
    start_reverse, stop_reverse, genes_r] = get_bsub_genes();

% find genes in the roi. 
if strand_region
    ind_genes = find( (start_forward>min(x_rec)) & ...
        (stop_forward<max(x_rec)));
else
    min_gene_other = genome_size-max(x_rec);
    max_gene_other = genome_size-min(x_rec);

    ind_genes = find(start_reverse>min_gene_other & ...
        stop_reverse<max_gene_other);
end

if ~strand_operon
    starts = genome_size-end_gene_operon+1;
    stops = genome_size-start_gene_operon+1;
    
    starts_roi = genome_size-stop_reverse(ind_genes)+1;
    stops_roi = genome_size-start_reverse(ind_genes)+1;
    genes_roi = {genes_r{ind_genes}};
else
    starts = start_gene_operon;
    stops = end_gene_operon;
    
    starts_roi = start_forward(ind_genes);
    stops_roi = stop_forward(ind_genes);
    genes_roi = {genes_f{ind_genes}};
end




w = 10;
dw = 2;
l = 100;

subplot(9,1,7);
hold on
for i = 1:length(ind_genes)
    [x,y] = arrow_function(starts_roi(i),stops_roi(i),w,dw,l);
    patch(x,y,'w')
    text(mean([starts_roi(i) stops_roi(i)]),w/2,genes_roi{i},...
        'HorizontalAlignment','center','Color','k','FontSize',12,'FontAngle','italic')
end

set(gca,'XLim',[(min(x_rec)-buffer) (min(x_rec)+dx_max+buffer)],'YLim',[-0.5*w 1.1*w+dw])
axis off



%% isoform stuff


% color scale for isoforms
isoform_levels_norm = isoform_levels/max(max(isoform_levels));
isoform_levels_norm(isoform_levels_norm<0) = 0;

min_level = log10(min(min(isoform_levels_norm(isoform_levels_norm>0))));
min_color = max([min_level min_color]);

n_colors = 500;
color_index = zeros(size(isoform_levels));
for i = 1:length(x_5)
    for j = 1:length(x_3)
        if log10(isoform_levels_norm(i,j))<min_color
            color_index(i,j) = 1;
        elseif log10(isoform_levels_norm(i,j))>max_color
            color_index(i,j)=n_colors;
        else
            color_index(i,j) = ceil(n_colors*(log10(isoform_levels_norm(i,j))-min_color)/(max_color-min_color))+1;
        end
    end
end
c = flipud(gray(n_colors+2));



% plotting isoforms
subplot(9,1,8:9);
hold on
counter = 1;
width = 6.5;
for i = 1:length(x_5)
    for j = 1:length(x_3)
        if log10(isoform_levels_norm(i,j))>min_color
            
            [x,y] = rectangle_points(x_5(i),x_3(j),-counter*10,width);
            patch(x,y,c(color_index(i,j),:),'LineWidth',0.5)
            
            plot([x_5(i) x_3(j)],-counter*10*[1 1],'-','LineWidth',6,'Color',c(color_index(i,j),:))
            text(x_3(j)+25,-counter*10,sprintf('%.2f',isoform_levels_norm(i,j)),'FontSize',8,...
                'VerticalAlignment','baseline')

            counter=counter+1;
        end
    end
end
colormap(c)
caxis([min_color max_color])
set(gca,'XLim',[(min(x_rec)-buffer) (min(x_rec)+dx_max+buffer)])
set(gca,'YLim',[-counter*10 0])
axis off
