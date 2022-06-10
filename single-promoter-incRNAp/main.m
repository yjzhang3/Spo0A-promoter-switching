%% start with just one vmax

% Ps
n_strain = 8;
ind = 1:8;
ii = 's';
vmax_array = struct();
group_array.g1 = 1:8;
file = 'May27-Psonly-1vmax.txt';

[bounds,final_energyi,vmax_array,diff] = estimate_energy(ii,n_strain,ind,file,group_array);
save('May27-Ps.mat')

% Pv
clear all
n_strain = 8;
ind = 1:8;
ii = 'v';
vmax_array = struct();
group_array.g1 = 1:8;
file = 'May27-Pvonly-1vmax.txt';

[bounds,final_energyi,vmax_array,diff] = estimate_energy(ii,n_strain,ind,file,group_array);
save('May27-Pv.mat')

%% plug in the vmax from 123*, and treat the rest with one vmax
% Ps
n_strain = 8;
ind = 1:8; % this is the index to choose dataset
ii = 's';
group_array.g1 = 1;
group_array.g2 = 2:8; % this is the index wrt bins
vmax_array.g1 = 762.8649;
file = 'May27-Psonly-2vmax.txt';

[bounds,final_energyi,vmax_array,diff] = estimate_energy(ii,n_strain,ind,file,group_array,vmax_array);
save('May27-Ps-2vmax.mat')

% Pv
clear all
n_strain = 8;
ind = 1:8;
ii = 'v';
group_array.g1 = 1;
group_array.g2 = 2:8; % this is the index wrt bins
vmax_array.g1 = 352.9209;
file = 'May27-Pvonly-2vmax.txt';

[bounds,final_energyi,vmax_array,diff] = estimate_energy(ii,n_strain,ind,file,group_array,vmax_array);
save('May27-Pv-2vmax.mat')

%% plug in the vmax from 123*, and treat the rest according to the results above
% worse fits will be put into one group
% but we still globally fit all 8 strains

% Ps
% WT, 1*(0xx),3*(xx0),23*(x00) fit worse
% hey, but a mutation type doesn't correspond to ONE configuration. For
% examplee, 1* can be the sum of 111 and 011! And WT is the sum of all 8 of
% configurations...

n_strain = 8;
ind = 1:8; % this is the index to choose dataset
ii = 's';
group_array.g1 = 1; % 000
group_array.g2 = [2,3,4]; % 001,010,011
group_array.g3 = [5,6,7,8]; % 111,110,101,100

vmax_array.g1 = 762.8649;
% vmax_array = struct();
file = 'May29-Psonly-3vmax.txt';

[bounds,final_energyi,vmax_array,diff] = estimate_energy(ii,n_strain,ind,file,group_array,vmax_array);
save('May29-Psonly-3vmax.mat')

% Pv
% WT(111),2*(x0x),12*(00x),3*(xx0),23*(x00) fit worse
clear all
n_strain = 8;
ind = 1:8;
ii = 'v';
group_array.g1 = 1;
group_array.g2 = [2,3,4]; % 001,010,011
group_array.g3 = [5,6,7,8]; % 100,101,110,111

vmax_array.g1 = 352.9209; 
% vmax_array = struct();
file = 'May29-Pvonly-3vmax.txt';

[bounds,final_energyi,vmax_array,diff] = estimate_energy(ii,n_strain,ind,file,group_array,vmax_array);
save('May29-Pvonly-3vmax.mat')

%% plug in the vmax from 123*, and treat the rest according to the results above
% worse fits will be put into one group
% but we still globally fit all 8 strains

% Ps
% WT, 1*(0xx),3*(xx0),23*(x00) fit worse
% hey, but a mutation type doesn't correspond to ONE configuration. For
% examplee, 1* can be the sum of 111 and 011! And WT is the sum of all 8 of
% configurations...

clear all
n_strain = 8;
ind = 1:8; % this is the index to choose dataset
ii = 's';
group_array.g1 = 1:4; % 000,001,010,011
group_array.g2 = 5:8; % 100,101,110,111

vmax_array.g1 = 762.8649;
% vmax_array = struct();
file = 'May30-Psonly-2vmax.txt';

[bounds,final_energyi,vmax_array,diff] = estimate_energy(ii,n_strain,ind,file,group_array,vmax_array);
save('May30-Psonly-2vmax.mat')

% Pv
% WT(111),2*(x0x),12*(00x),3*(xx0),23*(x00) fit worse
clear all
n_strain = 8;
ind = 1:8;
ii = 'v';
group_array.g1 = 1:4;% 000,001,010,011
group_array.g2 = 5:8;  % 100,101,110,111

vmax_array.g1 = 352.9209; 
% vmax_array = struct();
file = 'May30-Pvonly-2vmax.txt';

[bounds,final_energyi,vmax_array,diff] = estimate_energy(ii,n_strain,ind,file,group_array,vmax_array);
save('May30-Pvonly-2vmax.mat')