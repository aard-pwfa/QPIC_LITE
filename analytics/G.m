SI_consts;

n_step = 1001;
ns = linspace(5e15,5e17,n_step);% plas densities

g = zeros(1,n_step);            % geo factor
sd = zeros(1,n_step);           % nom skin depth
wp = zeros(1,n_step);           % nom freq
ez = zeros(1,n_step);           % ez no geo
EZ = zeros(1,n_step);           % ez w/geo
chi = zeros(1,n_step);          % chi
w = zeros(1,n_step);            % w*chi
lam = zeros(1,n_step);          % hollow wavelength
sdh = zeros(1,n_step);          % hollow sd
xp = zeros(1,n_step);           % ez*exp
EZT = zeros(1,n_step);          % ez with beam


N = 2e10;
sz = 40e-6;
for i = 1:n_step
    
%n0 = 8e16;
n0=ns(i);
[omega_p, lambda_p, skin_depth, plasma_time, plasma_period, E0, beta_p] = plasma_parameters(n0);

ez(i) = SI_e*omega_p^2/(pi*SI_eps0*SI_c^2);

a = 130/skin_depth;
b = 170/skin_depth;

b0 = B0(a,b);
b3 = B3(a,b);

g(i) = b0/(a*(2*b3+a*b0));
EZ(i) = ez(i)*g(i);

sd(i) = skin_depth;
wp(i) = omega_p;
chi(i) = sqrt(2*b3/(2*b3+a*b0));
w(i) = omega_p*chi(i);
lam(i) = 2*pi*SI_c/(chi(i)*omega_p);
sdh(i) = (chi(i)*omega_p)/SI_c;
xp(i) = N*exp(-(chi(i)*wp(i)*sz/SI_c)^2/2);
EZT(i) = xp(i)*EZ(i);

end