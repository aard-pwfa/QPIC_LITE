function plot_obj = PLOT_BEAM(beam,dim1,dim2,proj,nbin,d1lim,d2lim,plot_units,file_number,param_struct,plot_dir,save_plot,save_ext,figure_num)

num_str = num2str(file_number,'%04d');

cmap = custom_cmap;

if proj ~= 0
    if dim1 < 4
        PlotDim1 = param_struct.plasma.SD*beam(:,dim1)+proj*beam(:,dim1+3)./beam(:,6);
    end
    if dim2 < 4
        PlotDim2 = param_struct.plasma.SD*beam(:,dim2)+proj*beam(:,dim2+3)./beam(:,6);
    end
else
    if dim1 < 4
        PlotDim1 = param_struct.plasma.SD*beam(:,dim1);
    end
    if dim2 < 4
        PlotDim2 = param_struct.plasma.SD*beam(:,dim2);
    end
    if dim1 == 6
        PlotDim1 = (beam(:,dim1)-mean(beam(:,dim1)))/mean(beam(:,dim1));
    end
    if dim2 == 6
        PlotDim2 = (beam(:,dim2)-mean(beam(:,dim2)))/mean(beam(:,dim2));
    end
end

if isempty(d1lim)
    d1_edge = linspace(min(PlotDim1),max(PlotDim1),nbin);
else
    d1_edge = linspace(d1lim(1),d1lim(2),nbin);
end
if isempty(d2lim)
    d2_edge = linspace(min(PlotDim2),max(PlotDim2),nbin);
else
    d2_edge = linspace(d2lim(1),d2lim(2),nbin);
end

histmat = hist2(PlotDim1,PlotDim2,d1_edge,d2_edge);

figure(figure_num);
pcolor(d1_edge,d2_edge,histmat'); shading flat; box off;
colormap(cmap.wbgyr);

switch dim1
    case 1
        xlabel('X [\mum]','fontsize',14);
    case 2
        xlabel('Y [\mum]','fontsize',14);
    case 3
        xlabel('Z [\mum]','fontsize',14);
    case 4
        xlabel('X" [rad]','fontsize',14);
    case 5
        xlabel('Y" [rad]','fontsize',14);
    case 6
        xlabel('\delta','fontsize',14);
end

switch dim2
    case 1
        ylabel('X [\mum]','fontsize',14);
    case 2
        ylabel('Y [\mum]','fontsize',14);
    case 3
        ylabel('Z [\mum]','fontsize',14);
    case 4
        ylabel('X" [rad]','fontsize',14);
    case 5
        ylabel('Y" [rad]','fontsize',14);
    case 6
        ylabel('\delta','fontsize',14);
end

set(gca,'fontsize',14);
set(gcf,'color','w');

if save_plot
    saveas(gca,[plot_dir 'PS_d' num2str(dim1) '_d' num2str(dim2) '_' num_str save_ext]);
end
    
plot_obj = gca;