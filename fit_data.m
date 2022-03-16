function x = fit_data(nbd,TF_conc_t,mut_mat,real_data,lb,ub)
% input: real data, time dependent TF concentration, number of binding
% sites. Energyi is the parameter to find for

%% parameters
rng default

fun = @(p) objective_function(nbd,p,TF_conc_t,mut_mat,real_data);

nvars = 11; 

%% particle swarm
x = particleswarm(fun,nvars,lb,ub);

end