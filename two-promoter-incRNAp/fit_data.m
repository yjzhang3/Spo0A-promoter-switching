function [pars,M] = fit_data(nbd,TF_conc_t,RNApH_conc_t,RNApA_conc_t,mut_mat,real_data,lb,ub,n_strain,file)
% input: real data, time dependent TF concentration, number of binding
% sites, time-dependent RNAP concentration, mutation matrix

% can also specify lower and upper bound

% Energyi and vmax is the parameter to find for

% can specify how many strains are grouped together to optimize

%% parameters
rng default
iter = 8;

fun = @(p) objective_function(nbd,p,TF_conc_t,RNApH_conc_t,RNApA_conc_t,mut_mat,real_data,n_strain);

nvars = 15+24; % if only 3 binding boxes are included 

% energy of each binding site and vmax of each configuration

%% particle swarm
pars_all = zeros(iter,nvars);
diff_all = zeros(iter,1);
for n = 1:iter
    [x,fval,~,~] = particleswarm(fun,nvars,lb,ub);
    pars_all(n,:) = x;
    diff_all(n) = fval;
    create_parameter_file(file, x, fval);
end

[M,I] = min(diff_all);
pars =  pars_all(I,:);
end