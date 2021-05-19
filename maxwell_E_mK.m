function ret = maxwell_E_mK(E, T)
kB = 1.38e-23*1e-3; %J mK^-1
E=E.*kB;
f_E_J = 2*(kB*T)^-1.5*sqrt(E/pi) .*exp(-E./kB./T);
ret = f_E_J.*kB;