kB=1.38e-23 *1e-3; % J mK^-1 ; kg m^2 s^-2 mK^-1
u0=1.66e-27; % kg
%% Generate 1D Maxwell dist.
n1=51;n2=21;T=0.50; % 50 uK
plot_on=0;
[E,f_E,dE]= generate_prob_dens_energy(T,n1,n2,plot_on); %T and E in mK 
%%
syms y(t) 
w=1; %1 µm waist
m=133*u0;  %kg
U_0_T=1.1;  % mK
U_0 = -U_0_T/m*kB *1e-12 *1e12; %um^2 us^-2
fn_hm=diff(y, 2) == U_0*(y/w^2);
fn_gauss=diff(y, 2) == U_0*(y/w^2)*exp(-2*(y/w)^2);
% [V] = odeToVectorField(fn_gauss); % Gaussian
[V] = odeToVectorField(fn_hm); % Harmonic
M = matlabFunction(V,'vars', {'t','Y'}); % t in us, y in um

%%

t_interval=[0 400]; % us
sig=zeros(size(f_E));
for i=1:1:numel(f_E)
% for i=1:2
    E_i=E(i); %mK
    f_E_i=f_E(i);
%     x=sqrt(-log((U_0_T-E_i*1e-3)/(U_0_T))/2)*w; %um for Gaussian potential
%     x=sqrt(E_i/(U_0_T*1e3)/2)*w;   % µm  for Harmonic potential
    v_max=sqrt(2*E_i*kB/m); % m s^-1 or um us^-1
    x0=0;v0=v_max;
    initial_val = [v0 ; x0 ];
    sol = ode45(M,t_interval,initial_val);
    if plot_on==1
        figure(1)
        h=fplot(@(x1)deval(sol,x1,2),t_interval);
        sig(i)=sqrt(var(h.YData));
    else
        yDat=deval(sol,1:t_interval(end),2);
        sig(i)=sqrt(var(yDat));
    end
    xlabel('time (\mus)');
    ylabel('displacement (\mum)');
    title(['E = ' num2str(E_i) ' mK']);
end
dis= sum(sig.*f_E.*dE);
disp(['Mean displacement at T = ' num2str(T*1e3) 'µK in a ' num2str(U_0_T) ' mK trap is ' num2str(dis) ' µm']);
