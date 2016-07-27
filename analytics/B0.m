function x = B0(a,b)

b1 = besselk(0,a);
b2 = besseli(0,b);
b3 = besselk(0,b);
b4 = besseli(0,a);

x = b1.*b2 - b3.*b4;