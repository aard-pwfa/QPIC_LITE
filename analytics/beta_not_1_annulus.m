SI_consts;
n0 = 1e17;
[omega_p, lambda_p, skin_depth, plasma_time, plasma_period, E0, beta_p] = plasma_parameters(n0);

gamma = 40000;
beta = sqrt(1-gamma.^-2);

omega = linspace(0,2e13,20001);

k_p0 = omega./(gamma.*beta*SI_c);
k_p1 = omega./(gamma.*beta*SI_c)+omega_p/SI_c;

a = 200e-6;
b = 300e-6;

I_00a = besseli(0,k_p0*a);
I_10a = besseli(1,k_p0*a);
I_01a = besseli(0,k_p1*a);
I_11a = besseli(1,k_p1*a);

I_00b = besseli(0,k_p0*b);
I_10b = besseli(1,k_p0*b);
I_01b = besseli(0,k_p1*b);
I_11b = besseli(1,k_p1*b);

K_00a = besselk(0,k_p0*a);
K_10a = besselk(1,k_p0*a);
K_01a = besselk(0,k_p1*a);
K_11a = besselk(1,k_p1*a);

K_00b = besselk(0,k_p0*b);
K_10b = besselk(1,k_p0*b);
K_01b = besselk(0,k_p1*b);
K_11b = besselk(1,k_p1*b);

B00 = K_01a.*I_01b-K_01b.*I_01a;
B11 = K_11a.*I_11b-K_11b.*I_11a;
B10 = K_11a.*I_01b-K_01b.*I_11a;
B01 = K_01a.*I_11b-K_11b.*I_01a;

vareps = 1-omega_p^2./(omega.^2);

dispersion = k_p0.^2.*vareps.^2.*I_00a.*K_00b.*B11+k_p0.*k_p1.*vareps.*(I_00a.*K_10b.*B10+I_10a.*K_00b.*B01)+k_p1.^2.*I_10a.*K_10b.*B00;
[~,b] = min(abs(dispersion));
omega(b)

figure(1);
plot(omega,dispersion,'b',omega(b),0,'r*');

%%

dispersion = k_p0.*vareps.*I_00a.*K_11a+k_p1.*I_10a.*K_01a;
[~,b] = min(abs(dispersion));
omega(b)
figure(2);
plot(omega,dispersion,'b',omega(b),0,'r*');