function x = fit_data(real_data,TF_conc_t,RNAp_conc,nbd)
% input: real data, time dependent TF concentration, number of binding
% sites. Energyi is the parameter to find for

%% parameters
rng default
lb = zeros(11,1)-1E2; 
lb(1) = 1; % lower bound for RNAp concentration is nonzero

ub = zeros(11,1)+1E+2; 

fun = @(p) objective_function(real_data,TF_conc_t,RNAp_conc,nbd,p);

nvars = 11; 

%% particle swarm
x = particleswarm(fun,nvars,lb,ub);

end