% QuickPIC ANALYSIS SCRIPT
% SPENCER M.F. GESSNER 8/28

%clear all;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ANALYSIS AND PLOTTING FLAGS %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Name of folders containing data
% names = {'test2/','off_10/','off_20/','off_30/','off_40/','off_50/'};
% valsX = [0 10 20 30 40 50];
% scan_name = 'off_x/';

% names = {'test2/','off_10y/','off_20y/','off_30y/','off_40y/','off_50y/'};
% valsY = [0 10 20 30 40 50];
% scan_name = 'off_y/';

names = {'off_50/','off_15deg/','off_30deg/','off_45deg/','off_60deg/','off_75deg/', 'off_50y/'};
valsD = [0 15 30 45 60 75 90];
scan_name = 'off_angle/';

% Set Data and plotting directories
data_root = '/Users/sgess/Desktop/sims/data/qp_tars/Aug/';
plot_root = '/Users/sgess/Desktop/plots/QP/2014/Aug/';
plot_dir = [plot_root scan_name];

if ~exist(plot_dir,'dir')
    mkdir(plot_dir);
end

% Plot units, natural or real
plot_units = 'real';

% Save data?
save_plot = 0;
save_ext  = '.pdf';

% Get plasma and beam parameters
%load([data_dir 'param.mat']);
input_struct.plasma.SD = 18.7882;
input_struct.plasma.density = 8e16;
input_struct.plasma.field = 27.2;

% What data do you want to use?
QEB_type = 'QEB-XZ';
QEP1_type = 'QEP1-XZ';
FEZ_type = 'FEZ-XZ';

% Do beam analysis?
BEAM = 1;

% analyze last step
file_number = 100;

for i = 1:numel(names)
    
    data_dir = [data_root names{i}];
    
    %%%%%%%%%%%%%
    % LOAD DATA %
    %%%%%%%%%%%%%
    
    % Load beam and plasma densities
    beam_rho = LOAD_DATA(data_dir, QEB_type, file_number);
    plas_rho = LOAD_DATA(data_dir, QEP1_type, file_number);
    
    % Load fields
    field_EZ = LOAD_DATA(data_dir, FEZ_type, file_number);
    
    % Load axis and timing info
    [x_axis, z_axis] = LOAD_AXIS(data_dir, QEB_type, file_number);
    [iter, dt, time] = LOAD_TIME(data_dir, QEB_type, file_number);
    
    % Load 6D beam phase space
    if BEAM
        beam = LOAD_BEAM(data_dir,file_number);
        fit_dataX = BEAM_ANA(beam,1,2e6,128,[-700 700],input_struct);
        fit_dataY = BEAM_ANA(beam,2,2e6,128,[-700 700],input_struct);

        sigdX(i) = fit_dataX.sig;
        sigdY(i) = fit_dataY.sig;
        cendX(i) = fit_dataX.cen;
        cendY(i) = fit_dataY.cen;

        %display(fit_dataX.cen);
        %display(fit_dataY.cen);
    end
    
    
    
    %%%%%%%%%%%%%
    % PLOT DATA %
    %%%%%%%%%%%%%
    
    %PLOT_FUN(QEB_type,2,beam_rho,z_axis,x_axis,plot_units,file_number,input_struct,plot_dir,save_plot,save_ext,1);
    %PLOT_FUN(QEP1_type,2,plas_rho,z_axis,x_axis,plot_units,file_number,input_struct,plot_dir,save_plot,save_ext,2);
    %PLOT_FUN(FEZ_type,2,field_EZ,z_axis,x_axis,plot_units,file_number,input_struct,plot_dir,save_plot,save_ext,3);
    %PLOT_FUN(FEZ_type,1,field_EZ,z_axis,x_axis,plot_units,file_number,input_struct,plot_dir,save_plot,save_ext,4);
    %if BEAM; PLOT_BEAM(beam,1,2,2e6,128,[-700 700],[-700 700],plot_units,file_number,input_struct,plot_dir,save_plot,save_ext,15); end;
    %if BEAM; PLOT_BEAM(beam,1,6,0,128,[],[],plot_units,file_number,input_struct,plot_dir,save_plot,save_ext,8); end;
    
    pause(0.2);
end

% figure(15);
% plot(vals,sigX,'bs--',vals,sigY,'rs--',vals,cenX,'bo--',vals,cenY,'ro--','linewidth',2);
% set(gca,'fontsize',14);
% xlabel('Angle Offset [deg]','fontsize',14);
% ylabel('[\mum]','fontsize',14);
% l = legend('\sigma_{x}','\sigma_{y}','<X>','<Y>','location','northwest');
% set(l,'fontsize',14);
% set(gcf,'color','w');