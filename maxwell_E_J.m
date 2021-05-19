function ret = maxwell_E_J(E, T)
kB = 1.38e-23*1e-3; %J mK^-1
% E=E;
ret = 2*(kB*T)^-1.5*sqrt(E/pi) .*exp(-E./kB./T);