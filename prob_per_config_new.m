function p = prob_per_config_new(nbd, config, energyi,mut,TF_conc)

Z = Z_per_config_new(config,energyi,TF_conc);
Zall = Z_all_config_new(nbd,energyi,mut,TF_conc);

p = Z/Zall;

end 

