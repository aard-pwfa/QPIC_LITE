SI_consts;
n0 = 1e17;
[omega_p, lambda_p, skin_depth, plasma_time, plasma_period, E0, beta_p] = plasma_parameters(n0);

%gamma = 40000;
gamma = 2;
%gamma = 1:0.1:10;
beta = sqrt(1-gamma.^-2);

omega = linspace(0,0.5*omega_p,1000);
%omega = 0.7e13;

k_p0 = omega./(gamma.*beta*SI_c);
k_p1 = omega./(gamma.*beta*SI_c)+omega_p/SI_c;

a = 200e-6;
%a = (10:10:500)*1e-6;

I_10 = besseli(1,k_p0*a);
K_01 = besselk(0,k_p1*a);
I_00 = besseli(0,k_p0*a);
K_11 = besselk(1,k_p1*a);

dispersion = omega_p*(1+(k_p1.*I_10.*K_01)./(k_p0.*I_00.*K_11)).^(-1/2);

figure(1);
subplot(1,2,1);
plot(omega,omega,'b',omega,dispersion,'r','linewidth',2);
set(gca,'fontsize',12);
xlabel('\omega [rad/s]','fontsize',14);
ylabel('Dispersion [rad/s]','fontsize',14);
title('\gamma = 2','fontsize',14);

gamma = 40000;
%gamma = 2;
%gamma = 1:0.1:10;
beta = sqrt(1-gamma.^-2);

omega = linspace(0,0.5*omega_p,1000);
%omega = 0.7e13;

k_p0 = omega./(gamma.*beta*SI_c);
k_p1 = omega./(gamma.*beta*SI_c)+omega_p/SI_c;

a = 200e-6;
%a = (10:10:500)*1e-6;

I_10 = besseli(1,k_p0*a);
K_01 = besselk(0,k_p1*a);
I_00 = besseli(0,k_p0*a);
K_11 = besselk(1,k_p1*a);

dispersion = omega_p*(1+(k_p1.*I_10.*K_01)./(k_p0.*I_00.*K_11)).^(-1/2);

subplot(1,2,2);
plot(omega,omega,'b',omega,dispersion,'r','linewidth',2);
set(gca,'fontsize',12);
xlabel('\omega [rad/s]','fontsize',14);
ylabel('Dispersion [rad/s]','fontsize',14);
title('\gamma = 40000','fontsize',14);

%plot(a,dispersion)
%plot(omega,dispersion)
%plot(gamma,dispersion)