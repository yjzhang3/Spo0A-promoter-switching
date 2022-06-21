% fit single-promoter data in two batches and see if there is a fundamental
% conflict in particular energy values to make the fit only partially good
% for some mutants but not all of them

% '1vmax' is not accurate for some outputs. I forgot to change the name of
% the txt and mat file. Some outputs obviously have more vmax involved

%% Ps only

% batch 1: those that fit really well during global fit

n_strain = 4;
ind = [3,5,7:8];
ii = 's';
vmax_array = struct();
group_array.g1 = 1:8;
file = 'Jun7-Psonly-1vmax-b1.txt';

[bounds,final_energyi,vmax_array,diff] = estimate_energy(ii,n_strain,ind,file,group_array,vmax_array);
save('Jun7-Psonly-1vmax-b1.mat')

% batch 2: those that don't fit well with 2vmax (one unknown and one
% pre-opt)

% WT, 1*, 3*, 13* don't fit well
% (but we will use the same vmax from last fit)

n_strain = 4;
ind = [1:2,4,6];
ii = 's';
group_array.g1 = 1:8;
file = 'Jun7-Psonly-1vmax-b2.txt';

[bounds,final_energyi,vmax_array,diff] = estimate_energy(ii,n_strain,ind,file,group_array,vmax_array);
save('Jun7-Psonly-1vmax-b2.mat')

%% Ps only

% batch 3

% WT, 1*, 3*, 13* don't fit well, but WT seems to be mess up with this
% batch again. So take out WT.
% (but we will use the same vmax from last fit)

load('Jun7-Psonly-1vmax-b1.mat', 'vmax_array')

n_strain = 3;
ind = [2,4,6];
ii = 's';
group_array.g1 = 1:8;
file = 'Jun7-Psonly-1vmax-b3.txt';

[bounds,final_energyi,vmax_array,diff] = estimate_energy(ii,n_strain,ind,file,group_array,vmax_array);
save('Jun7-Psonly-1vmax-b3.mat')

%% Ps only

% batch 4: 1* only

% batch 3 overall is still not good. Try fitting just to 1*, and plot how
% the optimized energy does on 3*

load('Jun7-Psonly-1vmax-b1.mat', 'vmax_array')

n_strain = 1;
ind = 2;
ii = 's';
group_array.g1 = 1:8;
file = 'Jun7-Psonly-1vmax-b4.txt';

[bounds,final_energyi,vmax_array,diff] = estimate_energy(ii,n_strain,ind,file,group_array,vmax_array);
save('Jun7-Psonly-1vmax-b4.mat')


% use opt energies to generate data for 3* and 13*
load('promoter_activity_single.mat')

load('H_RNAP.mat', 'H_conc_t')

real_data_Ps = Ps_promoter_activity_mean(:,2:9); % exclude the real WT, instead WT is Ps only, no mutations
real_data_Ps_std = Ps_promoter_activity_std(:,2:9);

TF_conc_t = get_aps(2:10);

mut_mat = [[0,1,1];[1,1,0];[0,1,0]];
ind = [2,4,6];
nbd = 4;

real_data = real_data_Ps;
real_data_std = real_data_Ps_std;

title_name = {'Ps','1*','2*','3*','12*','13*','23*','123*'};

figure();
for kk = 1:length(ind)
    subplot(1,3,kk)
    errorbar(TF_conc_t,real_data(:,ind(kk)),real_data_std(:,ind(kk)),'LineStyle','none','LineWidth',2)
    hold on
    
    TR = time_dep_TR_new_wSigma(nbd,energy_copy,TF_conc_t,H_conc_t,mut_mat(kk,:),vmax_array,group_array);
    plot(TF_conc_t,TR,'LineWidth',2)
    xlabel('TF concentration')
    ylabel('transcription rate')
    
    ylim([0 1000])
    title(string(title_name(ind(kk))))
end
sgtitle('Using pre-opt energy from 1*')
saveas(gcf,'Jun7-b4-sim.jpeg')

%% another hypothesis is vmax for 1*,3*, and 13* is unique and different from the rest that fit weell
% judging from batch 4, 1* by itself fits really poorly. So the vmax from
% 123* fit is probably not good enough. Try to use more vmax (why not
% abondon it? because we know 123* vmax works well for other mutation types
% from previous runs!!)

% batch 5: 1* only and find a different vmax other than that used for 123*. Rather, use a
% different vmax. Specifically, group 1 includes 000, 100, 101,110, 111(pre-optimized from
% 123*), group 2 includes 001,010,011 (configurations that only
% mutation 1 correspond to)

% then we will see if the energy optimized from 1* only can apply well to
% 3* and 1*3* (1* and 3* share the same configurations that are 000,010)
% and 1* and 1*3* share 000,001

n_strain = 1;
ind = 2;
ii = 's';
group_array.g1 = [1,5,7,8];
group_array.g2 = [2:4,6];
load('Jun7-Psonly-1vmax-b1.mat', 'vmax_array')
file = 'Jun7-Psonly-1vmax-b5.txt';

[bounds,final_energyi,vmax_array,diff] = estimate_energy(ii,n_strain,ind,file,group_array,vmax_array);
save('Jun7-Psonly-1vmax-b5.mat')

% the fit is perfect for 1* only! So maybe now try use the energy values to
% generate data for 3* and 13*

% use opt energies to generate data for 3* and 13*
load('promoter_activity_single.mat')

load('H_RNAP.mat', 'H_conc_t')

real_data_Ps = Ps_promoter_activity_mean(:,2:9); % exclude the real WT, instead WT is Ps only, no mutations
real_data_Ps_std = Ps_promoter_activity_std(:,2:9);

TF_conc_t = get_aps(2:10);

mut_mat = [[0,1,1];[1,1,0];[0,1,0]]; %1*, 3*, 13*
ind = [2,4,6];
nbd = 4;

real_data = real_data_Ps;
real_data_std = real_data_Ps_std;

title_name = {'Ps','1*','2*','3*','12*','13*','23*','123*'};

figure();
for kk = 1:length(ind)
    subplot(1,3,kk)
    errorbar(TF_conc_t,real_data(:,ind(kk)),real_data_std(:,ind(kk)),'LineStyle','none','LineWidth',2)
    hold on
    
    TR = time_dep_TR_new_wSigma(nbd,final_energyi,TF_conc_t,H_conc_t,mut_mat(kk,:),vmax_array,group_array);
    plot(TF_conc_t,TR,'LineWidth',2)
    xlabel('TF concentration')
    ylabel('transcription rate')
    
    ylim([0 1000])
    title(string(title_name(ind(kk))))
end
sgtitle('Using pre-opt energy from 1*')
saveas(gcf,'Jun7-b5-sim.jpeg')

% generated data for 3* is really bad.

%% batch 6: I re-grouped the vmax (see notebook) and make sure the configuration unique to 1* and 3* are in their
% own group

% fit 1* and 3* together, using four vmax
n_strain = 2;
ind = [2,4];
ii = 's';
group_array.g1 = [1,6,8];
group_array.g2 = 3;
group_array.g3 = [2,4,5,7];
load('Jun7-Psonly-1vmax-b1.mat', 'vmax_array') %vmax_array.g1 is the pre-optimized one!
file = 'Jun7-Psonly-1vmax-b6.txt';

[bounds,final_energyi,vmax_array,diff] = estimate_energy(ii,n_strain,ind,file,group_array,vmax_array);
save('Jun7-Psonly-1vmax-b6.mat')

% both mutation types fit really well

%% batch 7: the above vmax grouping works well, so let's try fitting to 3 mutation types at the same time
% fit 1* and 3* together, using four vmax
n_strain = 3;
ind = [2,4,6];
ii = 's';
group_array.g1 = [1,6,8];
group_array.g2 = 3;
group_array.g3 = [2,4,5,7];
load('Jun7-Psonly-1vmax-b1.mat', 'vmax_array') %vmax_array.g1 is the pre-optimized one!
file = 'Jun7-Psonly-1vmax-b7.txt';

[bounds,final_energyi,vmax_array,diff] = estimate_energy(ii,n_strain,ind,file,group_array,vmax_array);
save('Jun7-Psonly-1vmax-b7.mat')

% all 3 works out well! What about we add more stuff in? all mutation types
% except for WT

%% batch 8: using the same vmax group above, but add all other mutation types except for WT
n_strain = 7;
ind = [2:8];
ii = 's';
group_array.g1 = [1,6,8];
group_array.g2 = 3;
group_array.g3 = [2,4,5,7];
load('Jun7-Psonly-1vmax-b1.mat', 'vmax_array') %vmax_array.g1 is the pre-optimized one!
file = 'Jun7-Psonly-1vmax-b8.txt';

[bounds,final_energyi,vmax_array,diff] = estimate_energy(ii,n_strain,ind,file,group_array,vmax_array);
save('Jun7-Psonly-1vmax-b8.mat')

%% batch 9: all mutation types and wild type

n_strain = 8;
ind = 1:8;
ii = 's';
group_array.g1 = [1,6,8];
group_array.g2 = 3;
group_array.g3 = [2,4,5,7];
load('H_RNAP', 'vmax') %vmax_array.g1 is the pre-optimized one!
vmax_array.g1 = vmax;
file = 'Jun7-Psonly-1vmax-b9.txt';

[bounds,final_energyi,vmax_array,diff] = estimate_energy(ii,n_strain,ind,file,group_array,vmax_array);
save('Jun7-Psonly-1vmax-b9.mat')

% all 3 works out well! What about we add more stuff in? all mutation types
% except for WT

