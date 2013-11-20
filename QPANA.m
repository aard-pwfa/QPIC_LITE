% QuickPIC ANALYSIS SCRIPT
% SPENCER M.F. GESSNER 9/8

%clear all;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ANALYSIS AND PLOTTING FLAGS %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Name of folder containing data
name = 'lb_ionize/';

% Set Data and plotting directories
%data_dir = '/Users/sgess/Desktop/FACET/2014/qp_tars/2013/Oct/SJG_004/';
%plot_dir = '/Users/sgess/Desktop/FACET/2014/qp_plots/2013/Oct/SJG_004/';

data_root = '/Volumes/PWFA_GoFlex/sgess/2014/';
plot_root = '/Users/sgess/Desktop/plots/QP/2014/';
data_dir = [data_root name];
plot_dir = [plot_root name];

if ~exist(plot_dir,'dir')
    mkdir(plot_dir);
end

% Plot units, natural or real
plot_units = 'real';

% Save data?
save_plot = 0;
save_ext  = '.pdf';

% Get plasma and beam parameters
load([data_dir 'param.mat']);
%input_struct = 0;

% What data do you want to use?
QEB_type = 'QEB-XZ';
QEP1_type = 'QEP1-XZ';
FEZ_type = 'FEZ-XZ';

% Do beam analysis?
BEAM = 1;

% Choose time steps to run over
file_number = 20;

for file_number = 20:20:60

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
end



%%%%%%%%%%%%%
% PLOT DATA %
%%%%%%%%%%%%%

PLOT_FUN(QEB_type,2,beam_rho,z_axis,x_axis,plot_units,file_number,input_struct,plot_dir,save_plot,save_ext,1);
PLOT_FUN(QEP1_type,2,plas_rho,z_axis,x_axis,plot_units,file_number,input_struct,plot_dir,save_plot,save_ext,2);
PLOT_FUN(FEZ_type,2,field_EZ,z_axis,x_axis,plot_units,file_number,input_struct,plot_dir,save_plot,save_ext,3);
PLOT_FUN(FEZ_type,1,field_EZ,z_axis,x_axis,plot_units,file_number,input_struct,plot_dir,save_plot,save_ext,4);
pause(0.2);
end
%save('../../../COMPARE/QP_eShort2.mat','ZZ','field_EZ');