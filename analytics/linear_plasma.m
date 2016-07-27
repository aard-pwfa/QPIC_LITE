cmap = custom_cmap;
SI_consts;

% plasma params
n0 = 8e16;
[omega_p, lambda_p, skin_depth, plasma_time, plasma_period, E0, beta_p] = plasma_parameters(n0);
kp=1/skin_depth;

% beam params
sr = 20;
sz = 15;
Nb = 4e8;

% axes
pts_z = 1000;
pts_x = 1000;
pts_r = 500;
l_z   = 22;
l_x   = 20;
l_r   = 10;
zz_cos = linspace(0,l_z*sz,pts_z+1);
zz_rho = linspace(-l_z*sz/2,l_z*sz/2,pts_z+1);
zz_pad = linspace(-l_z*sz/2,3*l_z*sz/2,2*pts_z+1);
dz     = zz_rho(2)-zz_rho(1);

xx     = linspace(-l_x*sr/2,l_x*sr/2,pts_x+1);
dx     = xx(2)-xx(1);
rr     = linspace(0,l_r*sr,pts_r+1);
dr     = rr(2)-rr(1);

% beam dist
rho_z    = 1/(sqrt(2*pi)*sz)*exp(-zz_rho.^2/(2*sz^2));
rho_zpad = 1/(sqrt(2*pi)*sz)*exp(-zz_pad.^2/(2*sz^2));
rho_r    = 1/(2*pi*sr^2)*exp(-rr.^2/(2*sr^2));
rho_x    = 1/(sqrt(2*pi)*sr)*exp(-xx.^2/(2*sr^2));
rho_xz   = rho_x'*rho_z;
rho_rz   = rho_r'*rho_z;
%% Calculations

% density calc
n1  = -kp*conv(rho_z,sin(kp*zz_cos),'full');
Dn1 = rho_x'*n1;
Dn_norm = Nb*Dn1*1e8/n0; % 1e8 converts um^2 to cm^2
Dn_OnAx = Dn_norm(pts_x/2+1,:);

% Ez calc
Wz = -(kp*1e6)^2*SI_e/SI_eps0*cos(kp*zz_cos);
Ez = dz*conv(Nb*rho_z,Wz,'full');

% Fr calc
Wf = -(kp*1e6)^2*SI_e/SI_eps0*sin(kp*zz_cos);
Fr = dz*conv(Nb*rho_z,Wf,'full');

% Ez R calc
R = zeros(size(rr));
R_f = zeros(size(rr));

dR = diff(rho_r)/dr;
dR(end+1) = 0;

in = zeros(size(rr));
out = zeros(size(rr));
in_f = zeros(size(rr));
out_f = zeros(size(rr));
for i = 1:numel(rr)
    
    r_lim = rr(i);
    r_lo = rr<=r_lim;
    r_hi = rr>r_lim;
    
    in(i) = sum(rho_r(r_lo).*besseli(0,rr(r_lo)/skin_depth).*besselk(0,r_lim/skin_depth).*rr(r_lo)*dr);
    out(i)= sum(rho_r(r_hi).*besseli(0,r_lim/skin_depth).*besselk(0,rr(r_hi)/skin_depth).*rr(r_hi)*dr);
    %in_f(i) = sum(rho_r(r_lo).*besseli(1,rr(r_lo)/skin_depth).*besselk(1,r_lim/skin_depth).*rr(r_lo)*dr);
    %out_f(i)= sum(rho_r(r_hi).*besseli(1,r_lim/skin_depth).*besselk(1,rr(r_hi)/skin_depth).*rr(r_hi)*dr);
    in_f(i) = sum(dR(r_lo).*besseli(1,rr(r_lo)/skin_depth).*besselk(1,r_lim/skin_depth).*rr(r_lo)*dr);
    out_f(i)= sum(dR(r_hi).*besseli(1,r_lim/skin_depth).*besselk(1,rr(r_hi)/skin_depth).*rr(r_hi)*dr);
    if i==1
        R(i) = out(i);
        R_f(i) = out_f(i);
    else
        R(i) = in(i)+out(i);
        R_f(i) = in_f(i)+out_f(i);
    end
    %R(i) = sum(rho_r(r_lo).*besseli(0,r(r_lo)/skin_depth).*besselk(0,r_lim/skin_depth).*r(r_lo)*dr)+...
    %       sum(rho_r(r_hi).*besseli(0,r_lim/skin_depth).*besselk(0,r(r_hi)/skin_depth).*r(r_hi)*dr);
end

RR = [fliplr(R), R(2:end)];
EzTot = RR'*Ez;

RF = [-fliplr(R_f), R_f(2:end)];
FrTot = RF'*Fr;
%%
figure(1);
imagesc(fliplr(zz_pad),xx,-Dn_norm);
colormap(cmap.bwr);  caxis([-max(abs(Dn_norm(:))) max(abs(Dn_norm(:)))]);
axis([-100 320 -120 120]);
set(gca, 'xdir','reverse');
xlabel('\xi [\mum]');
ylabel('r [\mum]');
h = colorbar;
ylabel(h,'n_1/n_0');
set(gca,'fontSize',16);

hold on;
plot(fliplr(zz_pad),10000*Dn_OnAx,'k','linewidth',2);
plot(5000*rho_x-100,xx,'k','linewidth',2)
hold off;

figure(2);
imagesc(fliplr(zz_pad),xx,-EzTot/1e9);
colormap(cmap.bwr);  caxis([-max(abs(EzTot(:)/1e9)) max(abs(EzTot(:)/1e9))]);
axis([-100 320 -120 120]);
set(gca, 'xdir','reverse');
xlabel('\xi [\mum]');
ylabel('r [\mum]');
h = colorbar;
ylabel(h,'E_z [GeV]');
set(gca,'fontSize',16);

hold on;
plot(fliplr(zz_pad),Ez/1e9,'k','linewidth',2);
plot(1500*RR-100,xx,'k','linewidth',2)
hold off;

figure(3);
imagesc(fliplr(zz_pad),xx,-FrTot/1e9);
colormap(cmap.bwr);  caxis([-max(abs(FrTot(:)/1e9)) max(abs(FrTot(:)/1e9))]);
axis([-100 320 -120 120]);
set(gca, 'xdir','reverse');
xlabel('\xi [\mum]');
ylabel('r [\mum]');
h = colorbar;
ylabel(h,'F_r/e [MT/m]');
set(gca,'fontSize',16);

hold on;
plot(-FrTot(:,501)/2e5,xx,'k','linewidth',2)
hold off;