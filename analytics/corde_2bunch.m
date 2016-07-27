% QuickPIC ANALYSIS SCRIPT
% SPENCER M.F. GESSNER 9/8

%clear all;
SI_consts;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ANALYSIS AND PLOTTING FLAGS %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Name of folder containing data
name = 'qpic_647/';

year = '2016/';

file_number = 001;

% z_lo = -150; 
% z_hi =350;
% beam_mag = 2;
% e_mag = 1200;
% shift_ind = 6;


save_plot= 0;
cmap = custom_cmap;

data_root = ['/Users/sgess/Desktop/sims/data/qp_tars/' year];
plot_root = ['/Users/sgess/Desktop/plots/QP/' year];
data_dir = [data_root name];
plot_dir = [plot_root name];
load([data_dir 'param.mat']);

if ~exist(plot_dir,'dir')
    mkdir(plot_dir);
end

% Plot units, natural or real
plot_units = 'real';

% Get plasma and beam parameters

% qpic 643
N1 = par.beam.Num_Particle(1);
N2 = par.beam.Num_Particle(2);
charge1 = par.beam.Charge(1);
charge2 = par.beam.Charge(2);
sz1 = par.beam.sigma_z(1)*1e-6;
sz2 = par.beam.sigma_z(2)*1e-6;
delta = (par.beam.Z_center(2)-par.beam.Z_center(1))*1e-6;
skin_depth = par.plasma.SD;
E0 = par.plasma.field;
n0 = par.plasma.Plasma_Density;
r_in = par.plasma.hollow_channel_inner_radius;
r_out = par.plasma.hollow_channel_outer_radius;
beam_mag = 10;
cur_mag = 5;

% What data do you want to use?
QEB_type = 'QEB-XZ';
QEP1_type = 'QEP1-XZ';
FEZX_type = 'FEZ-XZ';
FEYX_type = 'FEX-XZ';
FBXX_type = 'FBY-XZ';
% Do beam analysis?
BEAM = 1;

trans_ind = 2^(par.sim.INDX-1)+1;

% Load fields
field_EZX = 1000*E0*LOAD_DATA(data_dir, FEZX_type, file_number);
field_EYX = 1000*E0*LOAD_DATA(data_dir, FEYX_type, file_number);
field_BXX = 1000*E0*LOAD_DATA(data_dir, FBXX_type, file_number);
field_R = field_EYX-field_BXX;

EZ_line = field_EZX(trans_ind,:);

% Load axis and timing info
[iter, dt, time] = LOAD_TIME(data_dir, QEB_type, file_number);
[x_axis, z_axis] = LOAD_AXIS(data_dir, QEB_type, file_number);
xx = linspace(skin_depth*x_axis(1),skin_depth*x_axis(2),size(field_EZX,1));
xx_um = xx/1e6;
dx_um = xx_um(2)-xx_um(1);
zz = linspace(skin_depth*z_axis(1),skin_depth*z_axis(2),size(field_EZX,2));
dz = (zz(2)-zz(1))/1e6;

% Load beam and plasma densities
beam_rho = LOAD_DATA(data_dir, QEB_type, file_number);
plas_rho = LOAD_DATA(data_dir, QEP1_type, file_number);
tot_rho  = beam_mag*beam_rho + plas_rho;

% get beam and on-axis current
beam1 = LOAD_BEAM(data_dir,file_number,1);
beam2 = LOAD_BEAM(data_dir,file_number,2);

%%
shift_ind = 6;
Current1 = charge1*getCurrent(beam1,N1,n0,zz);
Current2 = charge2*getCurrent(beam2,N2,n0,zz);
Current = beam_mag*(Current1+Current2);





% ana stuff
a = r_in/skin_depth;
b = r_out/skin_depth;
if shift_ind > 0
    %[EZ_out, rho_b, shift_ind] = CompareTheory(zz,EZ_line,n0,N1,charge1,sz1,a,b,shift_ind);
    [EZ_out, rho_b, min_ind] = CompareTheory2(zz,EZ_line,n0,N1,charge1,sz1,N2,charge2,sz2,delta,a,b,shift_ind);
else
    %[EZ_out, rho_b, shift_ind] = CompareTheory(zz,EZ_line,n0,N,charge,sz,a,b);
    [EZ_out, rho_b, min_ind] = CompareTheory2(zz,EZ_line,n0,N1,charge1,sz1,N2,charge2,sz2,delta,a,b);

end
rho_c = beam_mag*(charge1*SI_c*SI_e*rho_b*(N1+N2)/(1000*sum(abs(rho_b)*dz)));

figure(1);
subplot(2,1,1);
% if charge < 0
%     subplot(2,2,1);
% else 
%     subplot(2,2,2);
% end
imagesc(zz,xx,tot_rho); c=colorbar('northoutside');colormap(flipud(cmap.bwr)); caxis([-1 1]);
axis xy;
set(gca,'xdir','reverse');
xlabel('Z [\mum]');
ylabel('Y [\mum]');
ylabel(c,'n/n_0');
set(gca,'fontsize',18);

%figure(2);
subplot(2,1,2);
plot(zz,cur_mag*Current,'b--',zz,EZ_line,'b',zz,cur_mag*(-rho_c),'r--',zz,-EZ_out,'r','linewidth',3); axis tight;
if charge1 < 0 
    %legend('Beam Profile, QuickPIC','E_z Field, QuickPIC','Beam Profile, Gaussian','E_z Field, Calculated','location','southeast');
else
    legend('Beam Profile, QuickPIC','E_z Field, QuickPIC','Beam Profile, Gaussian','E_z Field, Calculated','location','southwest');
end

%plot(zz,EZ_line,'b',zz,EZ_out,'r','linewidth',3);
%axis([z_lo z_hi -e_mag e_mag]);
set(gca,'xdir','reverse');
xlabel('Z [\mum]');
ylabel('E_z [MV/m]');
set(gca,'fontsize',18);

figure(2);
%subplot(2,1,1);
imagesc(zz,xx,field_EZX); c=colorbar('northoutside');colormap(flipud(cmap.bwr)); caxis([-200 200]);
axis xy;
set(gca,'xdir','reverse');
xlabel('Z [\mum]');
ylabel('Y [\mum]');
ylabel(c,'E_z [MV/m]');
set(gca,'fontsize',18);
% subplot(2,1,2);
figure(3);
imagesc(zz,xx,field_R); c=colorbar('northoutside');colormap(flipud(cmap.bwr)); %caxis([-200 200]);
axis xy;
set(gca,'xdir','reverse');
xlabel('Z [\mum]');
ylabel('Y [\mum]');
ylabel(c,'E_z [MV/m]');
set(gca,'fontsize',18);