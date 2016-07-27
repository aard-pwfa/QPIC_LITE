n0 = 2.5e16;

[omega_p, lambda_p, skin_depth, plasma_time, plasma_period, E0, beta_p] = plasma_parameters(n0);
sz = 30E-6;

zz_cos = linspace(0,20*sz,1001);
zz_rho = linspace(-10*sz,10*sz,1001);
zz_pad = linspace(-10*sz,30*sz,2001);
dz = zz_rho(2)-zz_rho(1);

N = 7E9;

rho = N/(sqrt(2*pi)*sz)*exp(-zz_rho.^2/(2*sz^2));
rho_pad = N/(sqrt(2*pi)*sz)*exp(-zz_pad.^2/(2*sz^2));

a = 275/skin_depth;
b = 325/skin_depth;
[Ez,Ez0,w0,chi] = holo_wake(n0,a,b,zz_cos);

Wz = -dz*conv(rho,Ez,'full');

%figure(1);
%plot(1e6*zz,dz*rho,'b');

%figure(2);
%plot(1e6*zz,Ez,'r');

figure(14);
%hold on;
%i_start = 1;
%i_end = 13000;
plot(1e6*zz_pad,100*dz*rho_pad/max(dz*rho_pad),'r',1e6*zz_pad,Wz/1e6,'b');



%%
% r = linspace(-15*sz,15*sz,15001);
% 
% b1 = besselk(0,a);
% b1r = besselk(0,r/(skin_depth*1e-6));
% b2 = besseli(0,b);
% b3 = besselk(0,b);
% b4 = besseli(0,a);
% b4r = besseli(0,r/(skin_depth*1e-6));
% 
% b0 = b1.*b2 - b3.*b4;
% br = b1r.*b2-b3.*b4r;
% Ezr = br/b0;
% Ezr(abs(r)<200e-6) = 1;
% Ezr(abs(r)>300e-6) = 0;
%Ezr(abs(r)>200 & abs(r)<300) = ;
%figure(5);
%plot(1e6*r,Ezr);

%%
% E_tot = Ezr'*Wz;
% 
% i_start = 300;
% i_end = 13000;
% r_start = 1800;
% r_end = 13200;
% cmap = custom_cmap;
% figure(6);
% imagesc(1e6*zz_pad(i_start:i_end)-10,1e6*r(r_start:r_end),real(E_tot(r_start:r_end,i_start:i_end))/1e9); caxis([-0.35 0.35]);
% colormap(cmap.bwr);
% 
% colorbar;