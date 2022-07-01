function [bounds,final_energyi,vmax_raw,diff,pars] = estimate_energy(ii,n_strain,ind,file,group_array)
%% this file finds the energy for smaller set of data
%% load data and choose strains
% ii tells us Pv or Ps
% n_strain tell us the number of strains to fit with
% ind tells us which strain to fit

% energyi = 1,2,3,p,12,13,1p,23,2p,3p

load('promoter_activity_single.mat')

load('H_RNAP.mat', 'H_conc_t')
load('A_RNAP.mat','A_conc_t')

real_data_Ps = Ps_promoter_activity_mean(:,2:9); % exclude the real WT, instead WT is Ps only, no mutations
real_data_Ps_std = Ps_promoter_activity_std(:,2:9);
real_data_Pv = Pv_promoter_activity_mean(:,2:9);
real_data_Pv_std = Pv_promoter_activity_std(:,2:9);

TF_conc_t = get_aps(2:10);

if ii == 'v'
    real_data = real_data_Pv;
    real_data_std = real_data_Pv_std;
    title_name = {'Pv','1*4*','2*4*','3*4*','124*','134*','234*','1234*'};
    RNAp_conc_t = A_conc_t;
else
    real_data = real_data_Ps;
    real_data_std = real_data_Ps_std;
    title_name = {'Ps','1*','2*','3*','12*','13*','23*','123*'};
    RNAp_conc_t = H_conc_t;
end


%% assign parameters
nbd = 4;
mut_mat = [[1,1,1];[0,1,1];[1,0,1];[1,1,0];[0,0,1];[0,1,0];[1,0,0];[0,0,0]]; 

%% set up bounds for energy
n_vars = 10+numel(fieldnames(group_array));
lb_overall = zeros(n_vars,1)-30;
ub_overall = zeros(n_vars,1)+20;

% lb_overall(1:3) = 1; % binding affinity of each site > 1 (nonzero and greater than promoter affinity)
% lb_overall(7) = 5; % must have one repressor

lb_overall(4) = 0;
ub_overall(4) = 0; % keep promoter energy fixed because we know it's embedded in concentration terms

% lb_overall(6) = 0;
% ub_overall(6) = 0; % TF bound at 0A1 and 3 probably won't have interaction energy

% vmax should be positive
lb_overall(11:end) = 0;
ub_overall(11:end) = 1000;

bounds.lb = lb_overall;
bounds.ub = ub_overall;

%% index smaller dataset
    
mut_mat_new = mut_mat(ind,:); 

real_data_new = real_data(:,ind);

%% start to fit data

% group and vmax array here specify which *strains* share the same vmax
[pars,diff] = fit_data(nbd,TF_conc_t,RNAp_conc_t,mut_mat_new,real_data_new,lb_overall,ub_overall,n_strain,group_array,ind,file);

%%
final_energyi = pars(1:10);
vmax_raw = pars(11:end);

fn = fieldnames(group_array);
vmax_per_strain = zeros(8,1);
for k=1:numel(fn)
    curr_ind = group_array.(fn{k});
    curr_vmax = vmax_raw(k);
    vmax_per_strain(curr_ind) = curr_vmax;
end
vmax_per_strain_final = vmax_per_strain(ind);

%% now see how the newer parameters do
figure();
for kk = 1:length(ind)
    subplot(4,2,kk)
    
    errorbar(TF_conc_t,real_data(:,ind(kk)),real_data_std(:,ind(kk)),'LineStyle','none','LineWidth',2)
    hold on
    
    TR = time_dep_TR_new_wSigma(nbd,final_energyi,TF_conc_t,RNAp_conc_t,mut_mat(ind(kk),:),vmax_per_strain_final(kk));
    plot(TF_conc_t,TR,'LineWidth',2)
    xlabel('TF concentration')
    ylabel('transcription rate')
    
    ylim([0 1000])
    title(string(title_name(ind(kk))))
end
saveas(gcf,append(file(1:end-4),'.jpeg'))
end
