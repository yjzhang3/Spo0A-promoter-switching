function Ztot = Z_overall_null(nbd,energyi,TF_conc,mRNA_conc)

Zon = Z_on(nbd,energyi,TF_conc,mRNA_conc);
Zoff = Z_off(nbd,energyi,TF_conc,mRNA_conc);
Ztot = Zon+Zoff;

end