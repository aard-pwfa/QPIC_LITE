% QuickPIC ANALYSIS SCRIPT
% SPENCER M.F. GESSNER 9/8

clear all;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ANALYSIS AND PLOTTING FLAGS %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Name of folder containing data
name = '30um/';

% Set Data and plotting directories
%data_dir = '/Users/sgess/Desktop/FACET/2014/qp_tars/2013/Oct/SJG_004/';
%plot_dir = '/Users/sgess/Desktop/FACET/2014/qp_plots/2013/Oct/SJG_004/';

data_root = '/Users/sgess/Desktop/sims/data/qp_tars/Aug/';
plot_root = '/Users/sgess/Desktop/plots/QP/2014/Aug/';
data_dir = [data_root name];
plot_dir = [plot_root name];

if ~exist(plot_dir,'dir')
    mkdir(plot_dir);
end

% Plot units, natural or real
plot_units = 'real';

% Save data?
save_plot = 0;
save_ext  = '.eps';

% Get plasma and beam parameters
%load([data_dir 'param.mat']);
n0 = 8e16;
[omega_p, lambda_p, skin_depth, plasma_time, plasma_period, E0, beta_p] = plasma_parameters(n0);
input_struct.plasma.SD = skin_depth;
input_struct.plasma.density = n0;
input_struct.plasma.field = E0;

% What data do you want to use?
QEB_type = 'QEB-XZ';
%QEB_type = 'QEB-YZ';
%QEB_type = 'QEB-XY';
QEP1_type = 'QEP1-XZ';
FEZX_type = 'FEZ-XZ';
%FEZY_type = 'FEZ-YZ';

% Do beam analysis?
BEAM = 1;

%for file_number = 10:10:100
file_number = 100;


% LOAD DATA %


% Load beam and plasma densities
beam_rho = LOAD_DATA(data_dir, QEB_type, file_number);
plas_rho = LOAD_DATA(data_dir, QEP1_type, file_number);

% Load fields
field_EZX = LOAD_DATA(data_dir, FEZX_type, file_number);
%field_EZY = LOAD_DATA(data_dir, FEZY_type, file_number);

% Load axis and timing info
[x_axis, z_axis] = LOAD_AXIS(data_dir, QEB_type, file_number);
xx = linspace(skin_depth*x_axis(1),skin_depth*x_axis(2),size(field_EZX,1));
zz = linspace(skin_depth*z_axis(1),skin_depth*z_axis(2),size(field_EZX,1));

[iter, dt, time] = LOAD_TIME(data_dir, QEB_type, file_number);



% Load 6D beam phase space
if BEAM
    beam = LOAD_BEAM(data_dir,file_number);
end
beam_z = hist(skin_depth*beam(:,3),zz);
[n,b] = hist(30*randn(1,100000),512);
%plot(skin_depth*zb,beam_z/max(beam_z),b,n/max(n));

%plot(zz,beam_z/max(beam_z),zz,E0*field_EZX(256,:));

%%

%%%%%%%%%%%%%
% PLOT DATA %
%%%%%%%%%%%%%

PLOT_FUN(QEB_type,2,-beam_rho,z_axis,x_axis,plot_units,file_number,input_struct,plot_dir,save_plot,save_ext,1);
PLOT_FUN(QEB_type,1,-beam_rho,z_axis,x_axis,plot_units,file_number,input_struct,plot_dir,save_plot,save_ext,42);

PLOT_FUN(QEP1_type,2,plas_rho,z_axis,x_axis,plot_units,file_number,input_struct,plot_dir,save_plot,save_ext,2);

PLOT_FUN(FEZX_type,2,field_EZX,z_axis,x_axis,plot_units,file_number,input_struct,plot_dir,save_plot,save_ext,30);
PLOT_FUN(FEZX_type,1,field_EZX,z_axis,x_axis,plot_units,file_number,input_struct,plot_dir,save_plot,save_ext,4);

%PLOT_FUN(FEZY_type,2,field_EZY,z_axis,x_axis,plot_units,file_number,input_struct,plot_dir,save_plot,save_ext,5);
%PLOT_FUN(FEZY_type,1,field_EZY,z_axis,x_axis,plot_units,file_number,input_struct,plot_dir,save_plot,save_ext,6);

if BEAM
    PLOT_BEAM(beam,1,2,2e6,128,[-700 700],[-700 700],plot_units,file_number,input_struct,plot_dir,save_plot,save_ext,15); 
    PLOT_BEAM(beam,3,6,0,128,[],[],plot_units,file_number,input_struct,plot_dir,save_plot,save_ext,8); 
end

PLOT_FUN(QEB_type,2,-10*beam_rho+plas_rho,z_axis,x_axis,plot_units,file_number,input_struct,plot_dir,save_plot,save_ext,7);


pause(0.2);