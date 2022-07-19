%% in this file, vmax is unique to each strain. 

%% for Ps, mutation on 3 would change the vmax
clear all
ii = 's';
ind = 1:8;
n_strain = 8;
group_array.g1 = [1,2,3,5]; % unaffected by the mutation
group_array.g2 = [4,6,7,8]; % affected by the mutation
file = 'Jun27-Ps-twogroup.txt';
[bounds,final_energyi,vmax_raw,diff,pars] = estimate_energy(ii,ind,file,group_array);
save('Jun27-Ps-twogroup.mat')

%% for Ps, mutation on at least 2 and/or 3 would change the vmax
clear all
ii = 's';
ind = 1:8;
n_strain = 8;
group_array.g1 = [1,2]; % unaffected by the mutation
group_array.g2 = [3,4,5,6,7,8]; % affected by the mutation
file = 'Jun29-Ps-twogroup.txt';
[bounds,final_energyi,vmax_raw,diff,pars] = estimate_energy(ii,ind,file,group_array);
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
[bounds,final_energyi,vmax_raw,diff,pars]= estimate_energy(ii,ind,file,group_array);
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
[bounds,final_energyi,vmax_raw,diff,pars] = estimate_energy(ii,ind,file,group_array);
save('Jun30-Pv-2group.mat')

%% we can do a confidence analysis on best parameters
% now store parameters in csv
% and do 100 iterations

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
file = 'Jul14-Ps-4group.csv';
[bounds,final_energyi,vmax_raw,diff,pars]= estimate_energy(ii,ind,file,group_array);
save('Jul14-Ps-4group.mat')

%% read csv and analyze parameters
plot_pars('Jul14-Ps-4group.csv',0.1)

%% for Pv, mutation on 1 would change the vmax
% but each type of mutation results in a differnet vmax
% WT,1,2,3,12,13,23,123
clear all
ii = 'v';
ind = 1:8;
n_strain = 8;
group_array.g1 = [1,3,4,7]; % unaffected by the mutation
group_array.g2 = [2,5,6,8]; % affected by the mutation
file = 'July17-Pv-2group.csv';
[bounds,final_energyi,vmax_raw,diff,pars] = estimate_energy(ii,ind,file,group_array);
save('July17-Pv-2group.mat')

%% read csv and analyze parameters
plot_pars('July17-Pv-2group.csv',0.1)
