%% Pv only

% batch 1: fit to all, but separate 000 and 100 by using different vmax for
% them
n_strain = 8;
ind = 1:8;
ii = 'v';
vmax_array = struct();
group_array.g1 = [1:4,6:8];
group_array.g2 = 5;
file = 'Jun8-Pvonly-2vmax-b1.txt';

[bounds,final_energyi,vmax_array,diff] = estimate_energy(ii,n_strain,ind,file,group_array,vmax_array);
save('Jun8-Pvonly-2vmax-b1.mat')

%% Pv only

% batch 2: fit to 23* only and see if separating 000 and 100 by using different vmax for
% them works
n_strain = 1;
ind = 7;
ii = 'v';
vmax_array = struct();
group_array.g1 = [1:4,6:8];
group_array.g2 = 5;
file = 'Jun8-Pvonly-2vmax-b2.txt';

[bounds,final_energyi,vmax_array,diff] = estimate_energy(ii,n_strain,ind,file,group_array,vmax_array);
save('Jun8-Pvonly-2vmax-b2.mat')

% now the fit is better

%% batch 3: can we improve better? this is only fit one set of data! should be good
n_strain = 1;
ind = 7;
ii = 'v';
vmax_array = struct();
group_array.g1 = 1;
group_array.g2 = 5;
group_array.g3 = [2:4,6:8];
file = 'Jun8-Pvonly-3vmax-b3.txt';

[bounds,final_energyi,vmax_array,diff] = estimate_energy(ii,n_strain,ind,file,group_array,vmax_array);
save('Jun8-Pvonly-3vmax-b3.mat')

% now the fit is only slightly better. So use two groups is okay.

%% batch 4: let's add in mutation types that share the same configurations as 23*
n_strain = 2;
ind = [3,7]; % 2* and 23* together
ii = 'v';
vmax_array = struct();
group_array.g1 = [1:4,6:8];
group_array.g2 = 5;
file = 'Jun8-Pvonly-2vmax-b4.txt';

[bounds,final_energyi,vmax_array,diff] = estimate_energy(ii,n_strain,ind,file,group_array,vmax_array);
save('Jun8-Pvonly-2vmax-b4.mat')

%% batch 5: let's add in mutation types that share the same configurations as 23*
n_strain = 3;
ind = [3,4,7]; % 2*, 3*,and 23* together
ii = 'v';
vmax_array = struct();
group_array.g1 = [1:4,6:8];
group_array.g2 = 5;
file = 'Jun8-Pvonly-2vmax-b5.txt';

[bounds,final_energyi,vmax_array,diff] = estimate_energy(ii,n_strain,ind,file,group_array,vmax_array);
save('Jun8-Pvonly-2vmax-b5.mat')

%% batch 6: let's add in mutation types that share the same configurations as 23*
n_strain = 8;
ind = 1:8; % let's fit them altogether
ii = 'v';
vmax_array = struct();
group_array.g1 = [1:4,6:8];
group_array.g2 = 5;
file = 'Jun8-Pvonly-2vmax-b6.txt';

[bounds,final_energyi,vmax_array,diff] = estimate_energy(ii,n_strain,ind,file,group_array,vmax_array);
save('Jun8-Pvonly-2vmax-b6.mat')

