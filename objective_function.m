function diff = objective_function(real_data,TF_conc_t,RNAp_conc,nbd,p)

% energyi = zeros(10,1);
% energyi(1) = pars.a;
% energyi(2) = pars.b;
% energyi(3) = pars.c;
% energyi(4) = pars.p;
% energyi(5) = pars.ab;
% energyi(6) = pars.ac;
% energyi(7) = pars.ap;
% energyi(8) = pars.bc;
% energyi(9) = pars.bp;
% energyi(10) = pars.cp;

tspan = length(real_data)*3600;

% TF_conc_t = @(x) x^3+30*x;

sim = time_dep_TR(TF_conc_t,tspan,RNAp_conc,nbd,p);
ind = 1:3600:tspan;
sim_data = sim(ind); % pick end of every hour

% now sim_data should be the same dimension as real_data

diff = sum((real_data-sim_data).^2);

end