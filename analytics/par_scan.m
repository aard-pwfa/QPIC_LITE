n0 = linspace(1e16,1e17,100);

a = zeros(size(n0));
b = zeros(size(n0));
omega_p = zeros(size(n0));
lambda_p = zeros(size(n0));
skin_depth = zeros(size(n0));
Ez0 = zeros(size(n0));
w0 = zeros(size(n0));
chi = zeros(size(n0));


for i = 1:numel(n0)

[omega_p(i), lambda_p(i), skin_depth(i)] = plasma_parameters(n0(i));

a(i) = 275/skin_depth(i);
b(i) = 325/skin_depth(i);
[Ez0(i),w0(i),chi(i)] = holo_wake_pars(n0(i),a(i),b(i));

end