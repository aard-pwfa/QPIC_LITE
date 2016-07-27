n0 = 8e16;

[omega_p, lambda_p, skin_depth, plasma_time, plasma_period, E0, beta_p] = plasma_parameters(n0);
sz = 30E-6;

zz_cos = linspace(0,20*sz,1001);
zz_rho = linspace(-10*sz,10*sz,1001);
zz_pad = linspace(-10*sz,30*sz,2001);
dz = zz_rho(2)-zz_rho(1);

N = 5E9;

rho = N/(sqrt(2*pi)*sz)*exp(-zz_rho.^2/(2*sz^2));
rho_pad = N/(sqrt(2*pi)*sz)*exp(-zz_pad.^2/(2*sz^2));

a = 200/skin_depth;
b = 300/skin_depth;
Ez = holo_wake(n0,a,b,zz_cos);

Wz = -dz*conv(rho,Ez,'full');

%figure(1);
%plot(1e6*zz_pad,dz*rho_pad/max(dz*rho_pad),'r',1e6*zz_pad,Wz/1e9,'b');