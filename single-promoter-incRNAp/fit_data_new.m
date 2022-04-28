function [pars,M] = fit_data_new(nbd,TF_conc_t,RNAp_conc_t,mut_mat,real_data,lb,ub)
% input: real data, time dependent TF concentration, number of binding
% sites, time-dependent RNAP concentration, mutation matrix

% can also specify lower and upper bound

% Energyi is the parameter to find for

%% parameters
rng default
iter = 4;

fun = @(p) objective_function_new(nbd,p,TF_conc_t,RNAp_conc_t,mut_mat,real_data);

nvars = 11; % if only 3 binding boxes are included 

%% particle swarm
pars_all = zeros(iter,nvars);
diff_all = zeros(iter,1);
for n = 1:iter
    x = particleswarm(fun,nvars,lb,ub);
    pars_all(n,:) = x;
    diff_all(n) = fun(x);
end

[M,I] = min(diff_all);
pars =  pars_all(I,:);
end