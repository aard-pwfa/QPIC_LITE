
SI_consts;

Ns = zeros(1,7);
min1_qp = zeros(1,7);
max1_qp = zeros(1,7);
min2_qp = zeros(1,7);
min1_ana = zeros(1,7);
max1_ana = zeros(1,7);
min2_ana = zeros(1,7);
lam_qp = zeros(1,7);
lam_ana = zeros(1,7);

nums = {'1E09','2E09','4E09','8E09','1E10','15E10','2E10'};
div = [1 1 1 1 1 10 1];

save_plot = 1;


for i = 1:7


% Name of folder containing data
num = nums{i};
name = ['pN' num '/'];
N = str2num(num)/div(i);
Ns(i) = N;
scale = 1.8e14;
charge = 1;
cax = [-.16*N/1e9 .16*N/1e9];

data_root = '/Users/sgess/Desktop/sims/data/qp_tars/2015/';
plot_root = '/Users/sgess/Desktop/plots/QP/2015/May/';
data_dir = [data_root name];
plot_dir = [plot_root name];

if ~exist(plot_dir,'dir')
    mkdir(plot_dir);
end

% Plot units, natural or real
plot_units = 'real';


n0 = 7e16;
[omega_p, lambda_p, skin_depth, plasma_time, plasma_period, E0, beta_p] = plasma_parameters(n0);
input_struct.plasma.SD = skin_depth;
input_struct.plasma.density = n0;
input_struct.plasma.field = E0;

% What data do you want to use?
QEB_type = 'QEB-XZ';

QEP1_type = 'QEP1-XZ';
FEZX_type = 'FEZ-XZ';

% Do beam analysis?
BEAM = 1;

%for file_number = 10:10:100
file_number = 10;


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
xx_um = xx/1e6;
zz = linspace(skin_depth*z_axis(1),skin_depth*z_axis(2),size(field_EZX,2));
zz_um = zz/1e6;
dz = zz(2)-zz(1);
dz_um = dz/1e6;

[iter, dt, time] = LOAD_TIME(data_dir, QEB_type, file_number);

% ana stuff
sz = 40E-6;
rho = charge*N/(sqrt(2*pi)*sz)*exp(-zz_um.^2/(2*sz^2));

a = 130/skin_depth;
b = 170/skin_depth;
[Ez, Ez0, w0] = holo_wake(n0,a,b,zz_um);

Wz = dz_um*conv(fliplr(rho),Ez,'full');
Wz = fliplr(Wz);
lz = zz_um(end)-zz_um(1);
zz_um2 = lz+zz_um(2:end);
zz_pad = [zz_um zz_um2];

good_ind = zz_pad <250/1e6;
%good_ind = zz_pad <650/1e6;

r = xx/skin_depth;

b1 = besselk(0,r);
b2 = besseli(0,b);
b3 = besselk(0,b);
b4 = besseli(0,r);

x = b1.*b2 - b3.*b4;

b1 = besselk(0,a);
b2 = besseli(0,b);
b3 = besselk(0,b);
b4 = besseli(0,a);

b0 = b1.*b2 - b3.*b4;

Ezx = x/b0;

r_in = abs(r)<a;
r_out = abs(r)>b;
Ezx(r_in) = 1;
Ezx(r_out) = 0;
Ezz=-Wz(good_ind)/1e9;
zz_tot = 1e6*zz_pad(good_ind);
Ez_tot = (Ezx')*Ezz;

ez_qp = E0*field_EZX(256,:);
ez_ana = -Wz(good_ind)/1e9;

min1_qp(i) = min(ez_qp(zz<50));
max1_qp(i) = max(ez_qp(zz>50));
min2_qp(i) = min(ez_qp(zz>50));

min1_ana(i) = min(ez_ana(zz<50));
max1_ana(i) = max(ez_ana(zz>50));
min2_ana(i) = min(ez_ana(zz>50));

%Ns = [Ns N];
%qpm = [qpm max(abs(E0*field_EZX(256,:)))];
%spm = [spm max(abs(Wz(good_ind)/1e9))];

% stuff

figure(1);
plot(zz, 2*charge*sum(beam_rho(256,:),1),'b--',zz,ez_qp,'b',  zz,2*charge*rho/scale,'r--',1e6*zz_pad(good_ind)-5,ez_ana,'r','linewidth',2);
axis tight;
set(gca, 'xdir','reverse');
h= legend('Beam Profile, QuickPIC','E_z Field, QuickPIC','Beam Profile, Gaussian','E_z Field, Calculated','location','southeast');
set(h,'fontsize',15);
set(gca,'fontsize',14);
xlabel('Z [\mum]','fontsize',16);
ylabel('E_z [GV/m]','fontsize',16);
saveas(gcf,[plot_dir 'EZ1D.png']);
%plot(zz_pad, Wz/1e9);
%plot(zz,E0*field_EZX(256,:))



dfz = 1/(zz(end)-zz(1));
fz_max = 1/dz;
n_fft = 10001;
kz = linspace(0,fz_max,n_fft);
yqp = fft(ez_qp,n_fft);
yana = fft(ez_ana,n_fft);

%figure(5);
%loglog(kz(1:round(n_fft/2)),abs(yqp(1:round(n_fft/2))),kz(1:round(n_fft/2)),abs(yana(1:round(n_fft/2))));
[r,t]=max(abs(yqp(1:round(n_fft/2))));
lam_qp(i)=1/kz(t);
[r,t]=max(abs(yana(1:round(n_fft/2))));
lam_ana(i)=1/kz(t);







figure(2);
imagesc(zz,xx,E0*field_EZX); set(gca, 'xdir','reverse');
cmap = custom_cmap;
colormap(cmap.bwr);
caxis(cax);
set(gca,'fontsize',14);
xlabel('Z [\mum]','fontsize',16);
ylabel('X [\mum]','fontsize',16);
t = colorbar('peer',gca);
set(get(t,'ylabel'),'String','E_z [GV/m]','fontsize',16);
saveas(gcf,[plot_dir 'EZ2D.png']);

% figure(3);
% imagesc(zz_tot,xx,real(Ez_tot)); set(gca, 'xdir','reverse');
% cmap = custom_cmap;
% colormap(cmap.bwr);
% caxis(cax);
% set(gca,'fontsize',14);
% xlabel('Z [\mum]','fontsize',16);
% ylabel('X [\mum]','fontsize',16);
% t = colorbar('peer',gca);
% set(get(t,'ylabel'),'String','E_z [GV/m]','fontsize',16);


figure(4);
mp=max(abs(plas_rho(:)));
mb=max(abs(beam_rho(:)));
sb= mp/(2*mb);
imagesc(zz,xx,-sb*beam_rho-plas_rho); set(gca, 'xdir','reverse');
%imagesc(zz,xx,-1*plas_rho); set(gca, 'xdir','reverse');
cmap = custom_cmap;
colormap(cmap.bwr);
%cm1 = flipud(bone(64));
%cm2 = hot(64);
%d = colormap([cm2;cm1]);
%colormap hot;

% hold on;
% imagesc(zz,xx,beam_rho); set(gca, 'xdir','reverse');
% colormap gray;

caxis([-1.0 1.0]);
set(gca,'fontsize',14);
xlabel('Z [\mum]','fontsize',16);
ylabel('X [\mum]','fontsize',16);
%t = colorbar('peer',gca);
%set(get(t,'ylabel'),'String','n/n_0','fontsize',16);
saveas(gcf,[plot_dir 'RHO.png']);
end

%%
figure(5);
plot(Ns,lam_ana,'r--',Ns,lam_qp,'bo','linewidth',2); axis([1e9 2e10 0 500]);
h= legend('Theory','QuickPIC','location','southeast');
set(h,'fontsize',15);
set(gca,'fontsize',14);
xlabel('N','fontsize',16);
ylabel('\lambda [\mum]','fontsize',16);
saveas(gcf,'~/Desktop/Plambda_scan.png');

%%
figure(6);
plot(Ns,min1_ana,'r--',Ns,min1_qp,'bo','linewidth',2); %axis([1e9 2e10 0 300]);
h= legend('Theory','QuickPIC','location','southeast');
set(h,'fontsize',15);
set(gca,'fontsize',14);
xlabel('N','fontsize',16);
ylabel('E_z [GeV/m]','fontsize',16);
saveas(gcf,'~/Desktop/Pmin1_scan.png');

%%
figure(7);
plot(Ns,max1_ana,'r--',Ns,max1_qp,'bo','linewidth',2); %axis([1e9 2e10 0 300]);
h= legend('Theory','QuickPIC','location','northeast');
set(h,'fontsize',15);
set(gca,'fontsize',14);
xlabel('N','fontsize',16);
ylabel('E_z [GeV/m]','fontsize',16);
saveas(gcf,'~/Desktop/Pmax1_scan.png');

%%
figure(8);
plot(Ns,min2_ana,'r--',Ns,min2_qp,'bo','linewidth',2); %axis([1e9 2e10 0 300]);
h= legend('Theory','QuickPIC','location','northeast');
set(h,'fontsize',15);
set(gca,'fontsize',14);
xlabel('N','fontsize',16);
ylabel('E_z [GeV/m]','fontsize',16);
saveas(gcf,'~/Desktop/Pmin2_scan.png');