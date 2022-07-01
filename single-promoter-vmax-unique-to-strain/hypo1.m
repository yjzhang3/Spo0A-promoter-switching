%% in this file, vmax is unique to each strain. 

%% for Ps, mutation on 3 would change the vmax
clear all
ii = 's';
ind = 1:8;
n_strain = 8;
group_array.g1 = [1,2,3,5]; % unaffected by the mutation
group_array.g2 = [4,6,7,8]; % affected by the mutation
file = 'Jun27-Ps-twogroup.txt';
[bounds,final_energyi,vmax_raw,diff] = estimate_energy(ii,n_strain,ind,file,group_array);
save('Jun27-Ps-twogroup.mat')

%% for Ps, mutation on at least 2 and/or 3 would change the vmax
clear all
ii = 's';
ind = 1:8;
n_strain = 8;
group_array.g1 = [1,2]; % unaffected by the mutation
group_array.g2 = [3,4,5,6,7,8]; % affected by the mutation
file = 'Jun29-Ps-twogroup.txt';
[bounds,final_energyi,vmax_raw,diff] = estimate_energy(ii,n_strain,ind,file,group_array);
save('Jun29-Ps-twogroup.mat')

%% for Ps, mutation on at least 2 and/or 3 would change the vmax
% but each type of mutation results in a differnet vmax
% WT,1,2,3,12,13,23,123
clear all
ii = 's';
ind = 1:8;
n_strain = 8;
group_array.g1 = [1,2]; % unaffected by the mutation
group_array.g2 = [3,5]; % affected by the mutation
group_array.g3 = [4,6];
group_array.g4 = [7,8];
file = 'Jun30-Ps-4group.txt';
[bounds,final_energyi,vmax_raw,diff] = estimate_energy(ii,n_strain,ind,file,group_array);
save('Jun30-Ps-4group.mat')

%% for Pv, mutation on 1 would change the vmax
% but each type of mutation results in a differnet vmax
% WT,1,2,3,12,13,23,123
clear all
ii = 'v';
ind = 1:8;
n_strain = 8;
group_array.g1 = [1,3,4,7]; % unaffected by the mutation
group_array.g2 = [2,5,6,8]; % affected by the mutation
file = 'Jun30-Pv-2group.txt';
[bounds,final_energyi,vmax_raw,diff] = estimate_energy(ii,n_strain,ind,file,group_array);
save('Jun30-Pv-2group.mat')