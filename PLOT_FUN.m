function plot_obj = PLOT_FUN(data_type,data_dim,data,x_axis,y_axis,plot_units,file_number,param_struct,plot_dir,save_plot,save_ext,figure_num)

num_str = num2str(file_number,'%04d');



%%%%%%%%%%%
% Density %
%%%%%%%%%%%

if strncmp(data_type,'QE',2)
    
    if strcmp(data_type(3),'B')
        particle = 'Beam';
    elseif strcmp(data_type(3),'P')
        particle = 'Plasma';
    end
    
    rho = data;
    
    if strcmp(plot_units,'natural')
        
        x_label = 'c/\omega_p';
        y_label = 'c/\omega_p';
        c_label = 'n_0';
        
        % Create axes
        ZZ = linspace(x_axis(1),x_axis(2),size(rho,1));
        XX = linspace(y_axis(1),y_axis(2),size(rho,2));
        
    elseif strcmp(plot_units,'real')
        
        x_label = '\mu m';
        y_label = '\mu m';
        c_label = 'cm^{-3}';
        
        % Create axes
        ZZ = linspace(x_axis(1),x_axis(2),size(rho,1))*param_struct.plasma.SD;
        XX = linspace(y_axis(1),y_axis(2),size(rho,2))*param_struct.plasma.SD;
        rho = rho*param_struct.plasma.density;
        
    end
        
    figure(figure_num);
    imagesc(ZZ,XX,rho);
    xlabel(x_label,'fontsize',16);
    ylabel(y_label,'fontsize',16);
    %caxis([-7 0]);
    colorbar;
    t = colorbar('peer',gca);
    set(get(t,'ylabel'),'String',c_label,'fontsize',16);
    title([particle ' Density'],'fontsize',16);
    
    if save_plot
        saveas(gca,[plot_dir particle '_rho_' num_str save_ext]);
    end
    
end



%%%%%%%%%%
% Fields %
%%%%%%%%%%

if strncmp(data_type,'F',1)
    
    if strcmp(data_type(2:3),'BX')
        comp = 'B_x';
    elseif strcmp(data_type(2:3),'BY')
        comp = 'B_y';
    elseif strcmp(data_type(2:3),'BZ')
        comp = 'B_z';
    elseif strcmp(data_type(2:3),'EX')
        comp = 'E_x';
    elseif strcmp(data_type(2:3),'EY')
        comp = 'E_y';
    elseif strcmp(data_type(2:3),'EZ')
        comp = 'E_z';
    end
    
    field = data;
    
    if strcmp(plot_units,'natural')
        
        x_label = 'c/\omega_p';
        y_label = 'c/\omega_p';
        c_label = 'm c \omega_p / e';
        
        % Create axes
        ZZ = linspace(x_axis(1),x_axis(2),size(field,1));
        XX = linspace(y_axis(1),y_axis(2),size(field,2));
        
    elseif strcmp(plot_units,'real')
        
        x_label = '\mu m';
        y_label = '\mu m';
        c_label = 'GV/m';
        
        % Create axes
        ZZ = linspace(x_axis(1),x_axis(2),size(field,1))*param_struct.plasma.SD;
        XX = linspace(y_axis(1),y_axis(2),size(field,2))*param_struct.plasma.SD;
        field = field*param_struct.plasma.field;
        
    end
    
    figure(figure_num);
    
    if data_dim == 1
        plot(ZZ,field(size(field,2)/2,:));
        xlabel(x_label,'fontsize',16);
        ylabel(c_label,'fontsize',16);
        title(comp,'fontsize',16);
    elseif data_dim == 2
        imagesc(ZZ,XX,field);
        xlabel(x_label,'fontsize',16);
        ylabel(y_label,'fontsize',16);
        colorbar;
        t = colorbar('peer',gca);
        set(get(t,'ylabel'),'String',c_label,'fontsize',16);
        title(comp,'fontsize',16);
    end
    
    if save_plot
        saveas(gca,[plot_dir comp '_field_' num_str save_ext]);
    end
    
end

plot_obj = gca;

