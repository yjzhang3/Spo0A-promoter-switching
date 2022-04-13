function Z = Z_per_config_new(config, energyi,TF_conc)

G0 = stand_energy(config,energyi);

[num_arr,conc_arr] = gen_config_prof(config,TF_conc,RNAp_conc);

num_TF = num_arr(1);
num_RNAp = num_arr(end);

conc_TF = conc_arr(1);
conc_RNAp = conc_arr(end);


Z = exp(-G0)*(conc_TF^num_TF)*(conc_RNAp^num_RNAp);

end