kB=1.38e-23; % kg m^2 s^-2 K^-1
syms y(t) 
w=1; %1 µm waist
m=133*1.6e-27;  %kg
U_0_T=1.1*1e-3;  % K
U_0 = -U_0_T/m *kB *1e-12 *1e12; %um^2 us^-2
% m=1;w=1;U_0=-1;
[V] = odeToVectorField(diff(y, 2) == U_0*(y/w^2)*exp(-2*(y/w)^2));

% [V] = odeToVectorField(m*diff(y, 2) == -U_0*(y/w^2))
M = matlabFunction(V,'vars', {'t','Y'});
t_interval=[0 200]; % us
initial_val=[1.1 0]; % um
sol = ode45(M,t_interval,initial_val);

figure(1)
h=fplot(@(x)deval(sol,x,1),t_interval);
xlabel('time (\mus)');
ylabel('displacement (\mum)');
hold on

[pks,pktimes]=findpeaks(h.YData,h.XData);
Period = mean(diff(pktimes));
freq =2*pi/Period % MHz
plot(h.XData,initial_val(1)*cos(freq*h.XData))