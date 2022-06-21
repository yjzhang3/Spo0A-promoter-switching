clear all
inds = 1:8;
indv = 1:8;
n_strains = 8;
n_strainv = 8;
group_arrays.g1 = [1,6,8];
group_arrays.g2 = 3;
group_arrays.g3 = [2,4,5,7];
load('H_RNAP', 'vmax') %vmax_array.g1 is the pre-optimized one!
vmax_array.g1 = vmax;
vmax_arrays = vmax_array;

vmax_arrayv = struct();
group_arrayv.g1 = [1:4,6:8];
group_arrayv.g2 = 5;

filename = 'Jun12-PvPs-besteachsetup.txt';

[bounds,vmax_arrays,vmax_arrayv,energyi_s,energyi_v,diff] = estimate_energy_together(inds,indv,n_strains,n_strainv,group_arrays,group_arrayv,vmax_arrays,vmax_arrayv,filename);
save('Jun12-PvPs-besteachsetup.mat')