function [E, f_E,dE] = generate_prob_dens_energy (T,n1,n2,plot_on)
    if nargin < 4
        plot_on=0;
    elseif nargin <3
        plot_on=0;
        n2=floor(n1/3);
    elseif nargin <2
        plot_on=0;
        n1=51;
        n2=floor(n1/3);
    end
    E=linspace(eps,T*5,n1);
    E=[E(1:end-1) linspace(E(end),T*10,n2) ] ;
    dE1=diff(E(1:2));
    dE2=diff(E(end-1:end));
    f_E = maxwell_E_mK(E,T);
    if plot_on
        figure(99)
        plot(E,f_E)
        xlabel('Energy (mK)');
        ylabel('Probability density f(E)');
        title(['For temperature of ' num2str(T) ' mK']);
        display(['Integral of (f(E)*deltaE) from E = [' num2str(E([1 end])) ']mK = ' num2str(sum([f_E(1:n1-1)*dE1 f_E(end-n2:end)*dE2])) ]);
    end
    dE=[ones(1,n1-1)*dE1 ones(1,n2)*dE2]; 
end
