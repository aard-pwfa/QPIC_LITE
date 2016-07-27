function [EZ_out, rho_b, min_ind] = CompareTheory2(zz,EZ_line,n0,N1,charge1,sz1,N2,charge2,sz2,delta,a,b,shift_ind)

ax_ind = numel(zz);
[~,ind] = min(abs(zz));
zz_cos = zz(ind:end)/1e6;
zz_rho = zz/1e6;
dz = zz_cos(2)-zz_cos(1);

rho_b1 = N1/(sqrt(2*pi)*sz1)*exp(-zz_rho.^2/(2*sz1^2));
rho_b2 = N2/(sqrt(2*pi)*sz2)*exp(-(zz_rho-delta).^2/(2*sz2^2));
rho_b = -charge1*rho_b1-charge2*rho_b2;

Ez = holo_wake(n0,a,b,zz_cos);
Wz = -charge1*dz*conv(rho_b,Ez,'full')/1e6;

diff_ind = numel(Wz)-ax_ind;
diff_sum = zeros(1,diff_ind);
for i = 1:diff_ind
    
    ind1 = i;
    ind2 = ax_ind+(i-1);
    Wz_comp = Wz(ind1:ind2);
    
    diff_sum(i) = sum((Wz_comp-EZ_line).^2);
end

if nargin < 13
    [~,min_ind] = min(diff_sum);
else
    min_ind = shift_ind;
end

EZ_out = Wz(min_ind:(ax_ind+(min_ind-1)));

