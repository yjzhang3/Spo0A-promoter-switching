function diff = objective_function(real_data,TF_conc_t,nbd,p)

RNAp_conc = p(1);
energyi = p(2:end);

real_data = real_data./max(real_data);

tspan = length(real_data)*3600;

% TF_conc_t = @(x) x^3+30*x;

sim = time_dep_TR(TF_conc_t,tspan,RNAp_conc,nbd,energyi);
ind = 1:3600:tspan;
sim_data = sim(ind); % pick end of every hour

% now sim_data should be the same dimension as real_data

diff = 1/length(real_data)*sum((real_data-sim_data).^2);

end