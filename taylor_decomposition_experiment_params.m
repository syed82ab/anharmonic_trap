FOV=0.38;
width=0.6;
width2=0.38;
addpath("D:\syeda\Repositories\SLM-controller\CPSWF\simple optimisation");
addpath("D:\syeda\Repositories\SLM-controller\CPSWF\p");


[E_r_fin, r_pswf_fin] = finding_optimal_spot_fwhm_etr(width2, FOV);

figure(50)
plot(r_pswf_fin,E_r_fin)
x=[-r_pswf_fin(end:-1:2); r_pswf_fin];
I_x=[E_r_fin(end:-1:2).^2; E_r_fin.^2]./E_r_fin(1).^2;
figure(49)
plot(x,I_x)
xlim([-1,1]);
xlabel('x (\lambda)')
%%
hold on
x_short=linspace(-1,1,1001);
a0=1; a2=-30; a4=3;
y2=a0+a2*x_short.^2+a4*x_short.^4;
plot(x_short,y2)
ylim([0 1.2])


%%
fun = @(x,xdata)x(1)+x(2)*xdata.^2+x(3)*xdata.^4+x(4)*xdata.^6+x(5)*xdata.^8+x(6)*xdata.^10+x(7)*xdata.^12;
x0=[a0 a2 a4 1 0 0 0];
figure(51)
plot(x_short,fun(x0,x_short),x,I_x)
xlim([-1,1]);
ylim([-1 1.2]);
%%
x_lim=0.4;
xdata=linspace(-x_lim,x_lim,501);
ydata=interp1(x,I_x,xdata);
figure(1)
hold on
x_fit=lsqcurvefit(fun,x0,xdata,ydata);
plot(x_short,fun(x_fit,x_short),x,I_x)
xlim([-1,1]);
ylim([-1 1.2]);

%%
figure(51)
plot(x,I_x,xdata,ydata,x_short,fun(x_fit,x_short))
xlim([-1,1])
ylim([0 2.2])
disp(['Polynomial coefficients for CPSWF are ' num2str(x_fit./x_fit(2))  ' at an x range of x =' char(0177) num2str(x_lim) ' ' char(0955)]);

%% Gaussian field

x_lim=0.9;
sig=1;
xdata=linspace(-x_lim,x_lim,501);
ydata=exp(-2*(xdata/sig).^2);
y_short=exp(-2*(x_short/sig).^2);
figure(52)
x0=[ 1 2 0 0 0 0 0];
xg_fit=lsqcurvefit(fun,x0,xdata,ydata);
plot(x_short,y_short,'b',xdata,ydata,'r',x_short,fun(xg_fit,x_short),'k')
ylim([0 1.1])
xg_fit./xg_fit(2);

disp(['Polynomial coefficients for Gaussian are ' num2str(xg_fit./xg_fit(2))  ' at an x range of x =' char(0177) num2str(x_lim) ' ' char(0955)]);
%%

figure
plot(x_short,y_short,x,E_x./E_r_fin(1).^2,x_short,1-x_short.^2)
xlim([-1 1])
xlabel('x (\mum)')
legend('Gaussian','CPSWF','Quadratic')