clearvars;
kB=1.38e-23 *1e-3; % J mK^-1 ; kg m^2 s^-2 mK^-1
u0=1.66e-27; % kg
%% Generate 1D Maxwell dist.
n1=51;n2=21;T=0.0010; % 
plot_on=0;disp_on=0;
[E,f_E,dE]= generate_prob_dens_energy(T,n1,n2,plot_on); %T and E in mK 
%%
syms y(t) omega epsil
w=1.; %1 µm gaussian waist
w_a=1.3; % 1.3 µm FWHM Airy disc waist
m=133*u0;  %kg
U_0_T=2.4;  % mK
U_0 = -U_0_T*kB *1e-12 *1e12; %um^2 us^-2
% k_hm=pi^2*U_0^2/2/w_airy^2;
% U_0_hm  = -U_0_T/m*kB *1e-12 *1e12; %um^2 us^-2
% U_0_gauss = -U_0_T/m*kB/w^2 *1e-12 *1e12; %um^2 us^-2
omega_r=2*pi*40e-3; % in MHz
k= m*omega_r^2;

fn_hm = diff(y, 2) == -k/m*(y)*(1+epsil*cos(omega*t));
% fn_gauss= diff(y, 2) == -U_0/m*gradient(exp(-2*(y/w)^2))*(1+epsil*cos(omega*t));
fn_gauss= diff(y, 2) == U_0/m*(4*(y)/(w^2))*exp(-2*(y/w)^2)*(1+epsil*cos(omega*t));
% fn_airy = diff(y, 2) == U_0/m*gradient(2*(besselj(1,pi*y/w_a)/(pi*y/w_a)).^2)*(1+epsil*cos(omega*t)); 
mk=pi/w_a;
fn_airy = diff(y, 2) == U_0./m.*(-8.*besselj(1,mk.*y)./((mk.*y).^2)./y.*(-2.*besselj(1,mk.*y)+mk.*y.*besselj(0,mk.*y))).*(1+epsil.*cos(omega.*t));
if exist('check','var')
    xx=linspace(0,10,201);
    U_g=exp(-2*(xx/w).^2);
    F_g=-(4*(xx)/(w^2)).*exp(-2*(xx/w).^2);
    figure(103)
    plot(xx,-gradient(U_g),xx,F_g*(xx(1)-xx(2)));
    legend('Gradient', 'Function used')
    xlabel('x(\mum)');ylabel('F*m (\muN)');
    title('Checking Gaussian Force from potential');
    
    U_a=(2*besselj(1,pi*xx/w_a)./(pi*xx/w_a)).^2;
%     F_a=-(8*w_a./(pi*xx.^2) .*( (w_a./(pi*xx)).^2.*(1+ pi*xx/w_a).*besselj(1,pi*xx/w_a).^2 - besselj(0,pi*xx/w_a).*besselj(1,pi*xx/w_a) ));
    mk=pi/w_a;
    F_a = -8*besselj(1,mk*xx)./(mk.*xx).^2 ./xx.*(besselj(1,mk*xx).*(-2) +mk*xx.*besselj(0,mk*xx));
%     F_a = -8*besselj(1,m*xx)./(m.*xx).^2 ./xx.*(-2*besselj(1,m*xx) +m*xx.*besselj(0,m*xx));

    figure(103)
    plot(xx,-gradient(U_a),xx,F_a*(xx(1)-xx(2)));
    legend('Gradient', 'Function used')
    xlabel('x(\mum)');ylabel('F*m (\muN)');
    title('Checking Airy disc Force from potential');
    
end

% [V] = odeToVectorField(fn_gauss); % Gaussian
% fn_name='gauss';
% [V] = odeToVectorField(fn_hm); % Harmonic
% fn_name='hm';
[V] = odeToVectorField(fn_airy); % Airy
fn_name='airy';

M = matlabFunction(V,'vars', {'t','Y','omega','epsil'}); % t in us, y in um

%%
% freq_scan = 2*pi*logspace(-2,0,51); 
freq_scan = linspace(50, 300, 201)*1e-3*2*pi; %MHz 2pi* 10-250kHz
epsil=0.04*ones(size(freq_scan));
dis=zeros(size(freq_scan));
sig=zeros(numel(freq_scan),numel(f_E));
for j=1:numel(freq_scan)
% for j=1

M1 = @(t,Y) M(t,Y,freq_scan(j),epsil(j));

t_interval=[0 400]; % us

for i=1:1:numel(f_E)
% for i=1:2
    E_i=E(i); %mK
    f_E_i=f_E(i);
%     x=sqrt(-log((U_0_T-E_i*1e-3)/(U_0_T))/2)*w; %um for Gaussian potential
%     x=sqrt(E_i/(U_0_T*1e3)/2)*w;   % µm  for Harmonic potential
    v_max=sqrt(2*E_i*kB/m); % m s^-1 or um us^-1
    x0=0;v0=-v_max;
    initial_val = [v0 ; x0 ];
    sol = ode45(M1,t_interval,initial_val);
    if plot_on==1
        figure(100)
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

h1=semilogy(freq_scan./(2*pi).*1e3,dis);

% h1=plot(freq_scan./(2*pi).*1e3,dis);
xlabel('\omega/(2\pi) (kHz)');ylabel('Mean displacement')
% trap_freq=sqrt(-4*U_0/m/w)*1e3; %kHz
trap_freq= sqrt(-U_0/m*pi^2/2/w_a)*1e3;
grid on
title({['Trap frequency = ' num2str(trap_freq/2/pi,3) ' kHz,'] ; 
       ['Atom Temperature = ' num2str(T) 'mK, Trap depth = ' num2str(U_0_T) ' mK']})
hold on

resonances=(2*trap_freq./(2*pi))./[1:4] ; %kHz
for jj=1:numel(resonances)
    xline(resonances(jj),'k--');
end
hold off
fig_fold= './new_figures/';
fig_filename = [fn_name '_U_0_' num2str(U_0_T) 'mK_atom_T_' num2str(T) 'mK_epsilon_' num2str(epsil(1)) '_mean_dis_vs_mod_freq.fig'];
savefig(gcf,[fig_fold fig_filename]);