%% same vmax, but with different promoter affinity for mutated strain
ii = 's';
ind = 1:8;
file = 'July1-Ps-onevmax-twodelta.txt';
group_vmax_array.g1 = 1:8; % share the same vmax
group_delta_array.g1 = [1,2,3,5]; % w/o 3* mutation
group_delta_array.g2 = [4,6,7,8]; % w/ 3* mutation
[bounds,diff,pars] = estimate_energy(ii,ind,file,group_vmax_array,group_delta_array);
save('July1-Ps-onevmax-twodelta.mat')

%% both vmax and delta have two groups
ii = 's';
ind = 1:8;
file = 'July4-Ps-twovmax-twodelta.txt';
group_vmax_array.g1 = [1,2,3,5]; 
group_vmax_array.g2 = [4,6,7,8];
group_delta_array.g1 = [1,2,3,5]; % w/o 3* mutation
group_delta_array.g2 = [4,6,7,8]; % w/ 3* mutation 
[bounds,diff,pars] = estimate_energy(ii,ind,file,group_vmax_array,group_delta_array);
save('July4-Ps-twovmax-twodelta.mat')

%% give 1*, 2*, and 12* mutation its own delta
ii = 's';
ind = 1:8;
group_vmax_array.g1 = 1:8; % one vmax
group_delta_array.g1 = [1,2]; % no mutation on 2* nor 3*
group_delta_array.g2 = [3,5]; % w/ 2* mutation
group_delta_array.g3 = [4,6]; % w/ 3* mutation
group_delta_array.g4 = [7,8]; % w/ 2* and 3* mutation simultaneously
file = 'July5-Ps-1vmax-4delta.txt';
[bounds,diff,pars] = estimate_energy(ii,ind,file,group_vmax_array,group_delta_array);
save('July5-Ps-1vmax-4delta.mat')

%% give 1*, 2*, and 12* mutation its own delta, and add more vmax
ii = 's';
ind = 1:8;
group_vmax_array.g1 = [1,2]; % no mutation on 2* nor 3*
group_vmax_array.g2 = [3,5]; % w/ 2* mutation
group_vmax_array.g3 = [4,6]; % w/ 3* mutation
group_vmax_array.g4 = [7,8]; % w/ 2* and 3* mutation simultaneously
group_delta_array.g1 = [1,2]; % no mutation on 2* nor 3*
group_delta_array.g2 = [3,5]; % w/ 2* mutation
group_delta_array.g3 = [4,6]; % w/ 3* mutation
group_delta_array.g4 = [7,8]; % w/ 2* and 3* mutation simultaneously
file = 'July5-Ps-4vmax-4delta.txt';
[bounds,diff,pars] = estimate_energy(ii,ind,file,group_vmax_array,group_delta_array);
save('July5-Ps-4vmax-4delta.mat')

%% for Pv, mutation on 1 would change the vmax
% but each type of mutation results in a differnet vmax
% WT,1,2,3,12,13,23,123
clear all
ii = 'v';
ind = 1:8;
group_vmax_array.g1 = 1:8;
group_delta_array.g1 = [1,3,4,7]; % unaffected by the mutation
group_delta_array.g2 = [2,5,6,8]; % affected by the mutation
file = 'July6-Pv-1vmax-2delta.txt';
[bounds,diff,pars] = estimate_energy(ii,ind,file,group_vmax_array,group_delta_array);
save('July6-Pv-1vmax-2delta.mat')