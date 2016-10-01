% load SI_constants
SI_consts;

% load plasma parameters
n0 = 8e16;
[omega_p, lambda_p, skin_depth, plasma_time, plasma_period, E0, beta_p] = plasma_parameters(n0);
kp = omega_p/SI_c;

% load geometric parameters
r_in = 240E-6;
r_out = 290E-6;

R_in = kp*r_in;
R_out = kp*r_in;
K1=besselk(1,R_in);
K2=besselk(2,R_in);
kap1 = kp^2*(K1/(R_in*K2))*(1+R_in*K1/(4*K2))^(-1);

% beam parameters
gamma = 39824;
sig_z = 35e-6;
N = 5.34e9;
I = N*SI_e*SI_c/sig_z;

% A parameter for BBU
A = 2*N*SI_re*kap1*sig_z/(gamma*r_in^2);
a_arg = A^(1/4);

% bunch position parameter for BBU
b = 0.5;
b_arg = sqrt(b);

% growth length
Lg = 2^(-3/2)*(N*SI_re*kap1*sig_z/(gamma*r_in^2))^(-1/2);
Lg0 = 2^(-5/2)*(N*SI_re*kap1*sig_z/(gamma*r_in^2))^(-1/2);

% Channel length
Lc = 0.073;
s = linspace(0,Lc,1000);
ds = s(2)-s(1);
s_arg = sqrt(s);

% growth functions
j0 = besselj(0,2*a_arg*s_arg*b_arg);
i0 = besseli(0,2*a_arg*s_arg*b_arg);
j1 = besselj(1,2*a_arg*s_arg*b_arg);
i1 = besseli(1,2*a_arg*s_arg*b_arg);

j0_Lg = besselj(0,s_arg/sqrt(Lg));
i0_Lg = besseli(0,s_arg/sqrt(Lg));
j1_Lg = besselj(1,s_arg/sqrt(Lg));
i1_Lg = besseli(1,s_arg/sqrt(Lg));

y0p = 0.4e-6;
y0 = 2.5*y0p;
y0p = 0;

y_out = y0/2*(j0+i0)+y0p*s_arg.*(j1+i1)/(2*a_arg*b_arg);
y_Lg = y0/2*(j0_Lg+i0_Lg)+y0p*s_arg.*(j1_Lg+i1_Lg)*(sqrt(Lg));
y_dri = y0+y0p*s;
y_app = y0/2*(besseli(0,s_arg/sqrt(Lg))+besselj(0,s_arg/sqrt(Lg)));
y_app2 = y0*(besseli(0,s_arg/sqrt(Lg)));


plot(100*s,1e6*y_out,'b',100*s,1e6*y_Lg,'r--',100*s,1e6*y_app,'g--')

yp_out = diff(y_out)/ds;
1e6*(yp_out(end)-y0p)

%%



