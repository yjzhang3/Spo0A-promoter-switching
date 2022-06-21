ii = 's';
ind = 1;
n_strain = 1;
file = 'Psonly-WT-1vmax.txt';
group_array.g1 = [1:8];
vmax_array = struct();


[bounds,final_energyi,vmax_array,diff] = estimate_energy(ii,n_strain,ind,file,group_array,vmax_array);
save('Psonly-WT-1vmax.mat')