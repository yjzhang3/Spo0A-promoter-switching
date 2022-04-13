% energyi = [A,B,C,Pv,AB,AC,APv,BC,BPv,CPv,Ps,AB,AC,APs,BC,BPs,CPs]
energyi = [10,10,10,1.5,0,0,20,0,-5,-5,2.7,0,0,-8,0,13,2];
nbd = 4;
mut = [1,1,1];
TF_conc_fun = @(x)x^3+30*x;
for ii = 1:500
    TF_conc_t(ii) = TF_conc_fun(ii);
end

%%
[TR_overall_v,TR_overall_s] = time_dep_TR_new(nbd,energyi,TF_conc_t,mut);