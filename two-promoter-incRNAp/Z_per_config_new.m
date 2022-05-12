function Z = Z_per_config_new(config, energyi,TF_conc,RNApH_conc,RNApA_conc)

G0 = stand_energy(config,energyi);

[num_arr,conc_arr] = gen_config_prof(config,TF_conc,RNApH_conc,RNApA_conc);

num_TF = num_arr(1);
num_RNApH = num_arr(2);
num_RNApA = num_arr(3);

conc_TF = conc_arr(1);
conc_RNApH = conc_arr(2);
conc_RNApA = conc_arr(3);

Z = exp(-G0)*(conc_TF^num_TF)*(conc_RNApH^num_RNApH)*(conc_RNApA^num_RNApA);

end