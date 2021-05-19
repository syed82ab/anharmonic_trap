kB=1.38e-23 *1e-3; % J mK^-1 ; kg m^2 s^-2 mK^-1
u0=1.66e-27; % kg
%% Generate 1D Maxwell dist.
n1=51;n2=21;T=0.50; % 50 uK
plot_on=0;disp_on;
[E,f_E,dE]= generate_prob_dens_energy(T,n1,n2,plot_on); %T and E in mK 
%%
syms y(t) omega epsil
w=1; %1 µm waist
m=133*u0;  %kg
U_0_T=1.1;  % mK
U_0 = -U_0_T/m*kB *1e-12 *1e12; %um^2 us^-2

fn_hm=diff(y, 2) == U_0*(y)*(1+epsil*cos(omega*t));
fn_gauss=diff(y, 2) == U_0*(y/w^2)*exp(-2*(y/w)^2)*(1+epsil*cos(omega*t));
% [V] = odeToVectorField(fn_gauss); % Gaussian
[V] = odeToVectorField(fn_hm); % Harmonic
M = matlabFunction(V,'vars', {'t','Y','omega','epsil'}); % t in us, y in um

%%
% freq_scan = 2*pi*logspace(-2,0,51); 
freq_scan = linspace(100, 550, 201)*1e-3; %MHz 100-550kHz
epsil=0.4*ones(size(freq_scan));
dis=zeros(size(freq_scan));
sig=zeros(numel(freq_scan),numel(f_E));
for j=1:numel(freq_scan)

M1 = @(t,Y) M(t,Y,freq_scan(j),epsil(j));

t_interval=[0 400]; % us

for i=1:1:numel(f_E)
% for i=1:2
    E_i=E(i); %mK
    f_E_i=f_E(i);
%     x=sqrt(-log((U_0_T-E_i*1e-3)/(U_0_T))/2)*w; %um for Gaussian potential
%     x=sqrt(E_i/(U_0_T*1e3)/2)*w;   % µm  for Harmonic potential
    v_max=sqrt(2*E_i*kB/m); % m s^-1 or um us^-1
    x0=0;v0=v_max;
    initial_val = [v0 ; x0 ];
    sol = ode45(M1,t_interval,initial_val);
    if plot_on==1
        figure(1)
        h=fplot(@(x1)deval(sol,x1,2),t_interval);
        sig(j,i)=sqrt(var(h.YData));
        xlabel('time (\mus)');
        ylabel('displacement (\mum)');
        title(['E = ' num2str(E_i) ' mK']);
    else
        yDat=deval(sol,1:t_interval(end),2);
        sig(j,i)=sqrt(var(yDat));
    end
end
dis(j)= sum(sig(j,i).*f_E.*dE);
if disp_on ==1
    disp(['Mean displacement at T = ' num2str(T*1e3) 'µK in a ' num2str(U_0_T) ' mK trap modulating at omega = ' num2str(freq_scan(j)) ' MHz is ' num2str(dis(j)) ' µm']);
end
end
%%
figure
h1=loglog(freq_scan./(2*pi).*1e3,dis);
xlabel('\omega/(2\pi) (kHz)');ylabel('Mean displacement')
trap_freq=sqrt(-U_0/w^2)*1e3; %kHz
title({['Trap frequency = ' num2str(trap_freq/2/pi,3) ' kHz,'] ; 
       ['Atom Temperature = ' num2str(T) 'mK, Trap depth = ' num2str(U_0_T) ' mK']})
hold on
resonances=(2*trap_freq./(2*pi))./[1:4] ; %kHz
for k=1:numel(resonances)
    xline(resonances(k),'r--');
end
hold off
fig_fold= 'figures\';
fig_filename = ['hm_U_0_' num2str(U_0_T) '_atom_T_' num2str(T) 'mean_dis_vs_mod_freq.fig'];
savefig(gcf,[fig_fold fig_filename]);