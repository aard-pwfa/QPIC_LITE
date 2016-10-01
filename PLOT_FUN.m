function [fig_num,plot_data,pAx1,pAx2] = PLOT_FUN(data_type,data,data_dim,data_dir,cmap,file_number,param_struct,plot_dir,save_plot,save_ext,figure_num)

num_str = num2str(file_number,'%04d');

if data_type(1)=='F'
    tag = data_type(2);
    vec = data_type(3);
    slice = data_type(5:6);
    if tag=='W'; data_type(2) = 'E'; end;
    [axis1, axis2] = LOAD_AXIS(data_dir, data_type, file_number);
    %data = LOAD_DATA(data_dir, data_type, file_number);
    if strcmp(slice,'XY'); data = data'; slice = 'YX'; end;
    if slice(1) ~= 'Z'
        axis1 = axis1 - mean(axis1);
    elseif slice(2) ~= 'Z'
        axis2 = axis2 - mean(axis2);
    end
    
    ax1 = linspace(axis1(1),axis1(2),size(data,1))*param_struct.plasma.SD;
    ax2 = linspace(axis2(1),axis2(2),size(data,2))*param_struct.plasma.SD;
    
    figure(figure_num);
    if data_dim == 1
        plot_data = param_struct.plasma.field*data(size(data,1)/2,:);
        pAx1 = ax2;
        plot(pAx1,plot_data,'k','linewidth',3); axis tight;
        if slice(2)=='Z'; set(gca,'xdir','reverse'); end;
        xlabel([slice(2) ' [\mum]']);
        if tag=='E'; ylabel([ tag '_{' vec '} [GV/m]']); end;
        if tag=='B'; ylabel([ tag '_{' vec '} [MT/m]']); end;
        if tag=='W'; ylabel([ tag '_{' vec '} [MT/m]']); end;
    elseif data_dim == 2
        plot_data = param_struct.plasma.field*data(:,size(data,2)/2);
        pAx1 = ax1;
        plot(pAx1,plot_data,'k','linewidth',3); axis tight;
        if slice(1)=='Z'; set(gca,'xdir','reverse'); end;
        xlabel([slice(1) ' [\mum]']);
        if tag=='E'; ylabel([ tag '_{' vec '} [GV/m]']); end;
        if tag=='B'; ylabel([ tag '_{' vec '} [MT/m]']); end;
        if tag=='W'; ylabel([ tag '_{' vec '} [MT/m]']); end;
    elseif data_dim == 3
        plot_data = param_struct.plasma.field*data;
        pAx1 = ax2;
        pAx2 = ax1;
        imagesc(pAx1,pAx2,plot_data); colormap(cmap);
        axis xy; 
        if slice(2) == 'Z'; set(gca,'xdir','reverse'); end;
        xlabel([slice(2) ' [\mum]']);
        ylabel([slice(1) ' [\mum]']);
        h=colorbar('northoutside'); 
        if tag=='E'; ylabel(h,[ tag '_{' vec '} [GV/m]']); end;
        if tag=='B'; ylabel(h,[ tag '_{' vec '} [MT/m]']); end;
        if tag=='W'; ylabel(h,[ tag '_{' vec '} [MT/m]']); end;
        val = max(abs(plot_data(:)));
        caxis([-val val]);
    end
    %title([ tag '_{' vec '}']);
    set(gca,'fontsize',24);
    
    if save_plot
        saveas(gca,[plot_dir tag '_' vec '_' slice num2str(data_dim) '_' num_str save_ext]);
    end
    
    
elseif data_type(1)=='Q'
    tag = data_type(3);
    slice = data_type((end-1):end);
    if tag=='T'; data_type(3) = 'B'; end;
    [axis1, axis2] = LOAD_AXIS(data_dir, data_type, file_number);
    %data = LOAD_DATA(data_dir, data_type, file_number);
    if strcmp(slice,'XY'); data = data'; slice = 'YX'; end;
    if slice(1) ~= 'Z'
        axis1 = axis1 - mean(axis1);
    elseif slice(2) ~= 'Z'
        axis2 = axis2 - mean(axis2);
    end
    
    ax1 = linspace(axis1(1),axis1(2),size(data,1))*param_struct.plasma.SD;
    ax2 = linspace(axis2(1),axis2(2),size(data,2))*param_struct.plasma.SD;
    figure(figure_num);
    if data_dim == 1
        plot_data = data(size(data,1)/2,:);
        pAx1 = ax2;
        plot(pAx1,plot_data,'k','linewidth',3); axis tight;
        if slice(2)=='Z'; set(gca,'xdir','reverse'); end;
        xlabel([slice(2) ' [\mum]']);
        ylabel(['n_{' tag '} [' num2str(param_struct.plasma.density,'%0.1E') ' cm^{-3}]']);
        
    elseif data_dim == 2
        plot_data = data(:,size(data,2)/2);
        pAx1 = ax1;
        plot(pAx1,plot_data,'k','linewidth',3); axis tight;
        if slice(1)=='Z'; set(gca,'xdir','reverse'); end;
        xlabel([slice(1) ' [\mum]']);
        ylabel(['n_{' tag '} [' num2str(param_struct.plasma.density,'%0.1E') ' cm^{-3}]']);

    elseif data_dim == 3
        plot_data = data;
        pAx1 = ax2;
        pAx2 = ax1;
        imagesc(pAx1,pAx2,plot_data); colormap(cmap);
        axis xy;
        if slice(2) == 'Z'; set(gca,'xdir','reverse'); end;
        xlabel([slice(2) ' [\mum]']);
        ylabel([slice(1) ' [\mum]']);
        h=colorbar('northoutside'); 
        ylabel(h,['n_{' tag '} [' num2str(param_struct.plasma.density,'%0.1E') ' cm^{-3}]']);
        caxis([-1 1]);
    end
    %if tag == 'B'; title('Beam Density'); end;
    %if tag == 'P'; title('Plasma Density'); end;
    %if tag == 'T'; title('Total Density'); end;
    set(gca,'fontsize',24);
    
    if save_plot
        saveas(gca,[plot_dir tag '_rho' '_' slice num2str(data_dim) '_' num_str save_ext]);
    end
    
end

fig_num = gca;