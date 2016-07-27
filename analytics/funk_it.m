cmap = custom_cmap;
SI_consts;
n0 = 8e16;

[omega_p, lambda_p, skin_depth, plasma_time, plasma_period, E0, beta_p] = plasma_parameters(n0);
kp = 1/skin_depth;

sz = 15;
%N = 0.5E10;
N = 4E8;
sr = 20;

zz_cos = linspace(0,20*skin_depth,1001);
zz_rho = linspace(-10*skin_depth,10*skin_depth,1001);
zz_pad = linspace(-10*skin_depth,30*skin_depth,2001);
xx     = linspace(-10*skin_depth,10*skin_depth,1001);
dz = zz_rho(2)-zz_rho(1);
dx = xx(2)-xx(1);

rho_z = 1/(sqrt(2*pi)*sz)*exp(-zz_rho.^2/(2*sz^2));
rho_zpad = 1/(sqrt(2*pi)*sz)*exp(-zz_pad.^2/(2*sz^2));

Wz = -(kp*1e6)^2*SI_e/SI_eps0*cos(kp*zz_cos);

Ez = dz*conv(rho_z,Wz,'full');

figure(1);
plot(zz_pad,Ez);

%%
r = linspace(0,10*skin_depth,501);
dr = r(2)-r(1);

rho_r = 1/(2*pi*sr^2)*exp(-r.^2/(2*sr^2));
R = zeros(size(r));
in = zeros(size(r));
out = zeros(size(r));
for i = 1:numel(r)
    
    r_lim = r(i);
    r_lo = r<=r_lim;
    r_hi = r>r_lim;
    
    in(i) = sum(rho_r(r_lo).*besseli(0,r(r_lo)/skin_depth).*besselk(0,r_lim/skin_depth).*r(r_lo)*dr);
    out(i)= sum(rho_r(r_hi).*besseli(0,r_lim/skin_depth).*besselk(0,r(r_hi)/skin_depth).*r(r_hi)*dr);
    if i==1
        R(i) = out(i);
    else
        R(i) = in(i)+out(i);
    end
    %R(i) = sum(rho_r(r_lo).*besseli(0,r(r_lo)/skin_depth).*besselk(0,r_lim/skin_depth).*r(r_lo)*dr)+...
    %       sum(rho_r(r_hi).*besseli(0,r_lim/skin_depth).*besselk(0,r(r_hi)/skin_depth).*r(r_hi)*dr);
end

figure(11);
plot(r,170*rho_r,'b',r,R,'k');
%%
figure(2);
rho_t = rho_r'*rho_z;
imagesc(rho_t);

Etot = N*R'*Ez;

Etot2 = flipud(Etot);

Eall = [Etot2; Etot];

figure(3);
imagesc(zz_pad,xx,Eall/1e9); colormap(cmap.bwr); colorbar;
set(gca,'xdir','reverse');
axis([-100 375 -80 80]);
caxis([-1 1]);

%rho_x = 1/(sqrt(2*pi)*sr)*exp(-xx.^2/(2*sr^2));
%rho_t = rho_r'*rho_z;
%rho_tpad = rho_x'*rho_zpad;
