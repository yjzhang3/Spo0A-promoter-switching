function Z = Z_perconfig(config,energyi,TF_conc,RNAp_conc)

Gtot = free_energy(config,energyi,TF_conc,RNAp_conc);

Z = exp(-Gtot);

end