%% axes

f=figure(1);
n1=2;n2=3;
pos=zeros(n1*n2,4);
for i=1:n1*n2
    subplot(n1,n2,i);
end
for i=1:n1*n2
    pos(i,:)=f.Children(i).Position;
end
for i=n1*n2:-1:1
    delete(f.Children(i));
end


fig_fold='./figures/';
fig1={'even_U_0_1.1mK_atom_T_0.005mK_epsilon_0.04_coeff_1_0_0_0_0_0_mean_dis_vs_mod_freq.fig';
      'even_U_0_1.1mK_atom_T_0.005mK_epsilon_0.04_coeff_1_-1_0.66_-0.32_0.11_-0.021_mean_dis_vs_mod_freq.fig';
      'even_U_0_1.1mK_atom_T_0.005mK_epsilon_0.04_coeff_1_-3.5_7_-8.1_5.3_-1.5_mean_dis_vs_mod_freq.fig' ;
      'even_U_0_1.1mK_atom_T_0.05mK_epsilon_0.04_coeff_1_0_0_0_0_0_mean_dis_vs_mod_freq.fig';
      'even_U_0_1.1mK_atom_T_0.05mK_epsilon_0.04_coeff_1_-1_0.66_-0.32_0.11_-0.021_mean_dis_vs_mod_freq.fig';
      'even_U_0_1.1mK_atom_T_0.05mK_epsilon_0.04_coeff_1_-3.5_7_-8.1_5.3_-1.5_mean_dis_vs_mod_freq.fig' ;
    };
for i=1:n1*n2
    f1(i)=openfig([fig_fold fig1{i}]);
end
%%
for i=1:n1*n2
    copyobj(f1(i).Children,f)
end
for i=1:n1*n2
    f.Children(i).Position = pos(i,:);
end