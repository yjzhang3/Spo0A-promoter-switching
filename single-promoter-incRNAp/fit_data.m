function [pars,M] = fit_data(nbd,TF_conc_t,RNAp_conc_t,mut_mat,real_data,lb,ub,n_strain)
% input: real data, time dependent TF concentration, number of binding
% sites, time-dependent RNAP concentration, mutation matrix

% can also specify lower and upper bound

% Energyi is the parameter to find for

%% parameters
rng default
iter = 4;

fun = @(p) objective_function(nbd,p,TF_conc_t,RNAp_conc_t,mut_mat,real_data,n_strain);

nvars = 10+8; % if only 3 binding boxes are included 

% energy of each binding site and vmax of each configuration

%% particle swarm
pars_all = zeros(iter,nvars);
diff_all = zeros(iter,1);
for n = 1:iter
    [x,fval,~,~] = particleswarm(fun,nvars,lb,ub);
    pars_all(n,:) = x;
    diff_all(n) = fval;
end

[M,I] = min(diff_all);
pars =  pars_all(I,:);
end