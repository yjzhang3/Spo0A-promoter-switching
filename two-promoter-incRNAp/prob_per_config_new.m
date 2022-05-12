function p = prob_per_config_new(nbd, config, energyi,mut,TF_conc,RNApH_conc,RNApA_conc)

Z = Z_per_config_new(config,energyi,TF_conc,RNApH_conc,RNApA_conc);
Zall = Z_all_config_new(nbd,energyi,mut,TF_conc,RNApH_conc,RNApA_conc);

p = Z/Zall;

end 

