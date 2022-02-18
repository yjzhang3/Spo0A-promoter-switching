function Z = Z_perconfig(config,energyi,TF_conc,mRNA_conc)

Gtot = free_energy(config,energyi,TF_conc,mRNA_conc);

Z = exp(-Gtot);

end