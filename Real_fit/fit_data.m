function [pars,diff] = fit_data(nbd,TF_conc_t,mut_mat,real_data,lb,ub)
% input: real data, time dependent TF concentration, number of binding
% sites. Energyi is the parameter to find for

%% parameters
rng default
iter = 12;

fun = @(p) objective_function(nbd,p,TF_conc_t,mut_mat,real_data);

nvars = nbd+nbd*(nbd-1)/2;

%% particle swarm
pars = zeros(iter,nvars);
diff = zeros(iter,1);
parfor n = 1:iter
    x = particleswarm(fun,nvars,lb,ub);
    pars(n,:) = x;
    diff(n) = fun(x);
end


end