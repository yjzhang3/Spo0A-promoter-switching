function [pars,M] = fit_data_together(nbd,TF_conc_t,H_conc_t,A_conc_t,mut_mats,mut_matv,real_datas,real_datav,lb,ub,n_strains,n_strainv)
% input: real data, time dependent TF concentration, number of binding
% sites, time-dependent RNAP concentration, mutation matrix

% can also specify lower and upper bound

% Energyi and vmax is the parameter to find for

% can specify how many strains are grouped together to optimize

%% parameters
rng default
iter = 4;

fun = @(p) objective_function_together(nbd,p,TF_conc_t,H_conc_t,A_conc_t,mut_mats,mut_matv,real_datas,real_datav,n_strains,n_strainv);
% nvars = 10+8; % if only 3 binding boxes are included 
nvars = 10+6+8+8;

% energy of each binding site and vmax of each configuration

%% particle swarm
pars_all = zeros(iter,nvars);
diff_all = zeros(iter,1);
for n = 1:iter
    [x,fval,~,~] = particleswarm(fun,nvars,lb,ub);
    pars_all(n,:) = x;
    diff_all(n) = fval;
    create_parameter_file('pars_together.txt', x, fval)
end

[M,I] = min(diff_all);
pars =  pars_all(I,:);
end