function [bounds,final_energyi,vmax_array,diff] = estimate_energy(ii,n_strain,ind,file,group_array,vmax_array)
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

fnv = fieldnames(vmax_array);

if isempty(fnv)
    n_vars = 10+numel(fieldnames(group_array));
end
if ~isempty(fnv)
    n_vars = 10+numel(fieldnames(group_array))-numel(fieldnames(vmax_array));
end

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
[pars,diff] = fit_data(nbd,TF_conc_t,RNAp_conc_t,mut_mat_new,real_data_new,lb_overall,ub_overall,n_strain,group_array,vmax_array,file);

%%
final_energyi = pars(1:10);
vmax = pars(11:end);

%%

fn = fieldnames(group_array);
ol = numel(fieldnames(vmax_array)); % original length of vmax array field
for k=1:numel(fn)
    curr_key = (fn{k});
    if ~isfield(vmax_array,curr_key)  % only if this group doesn't exist, then we have unknown vmax
        curr_value = vmax(k-ol);
        vmax_array.(fn{k}) = curr_value; % assign unknown parameters
    end
    % of the vmax's to be optimized the same as the order of the group
end

%% now see how the newer parameters do
figure();
for kk = 1:length(ind)
    subplot(4,2,kk)
    
    errorbar(TF_conc_t,real_data(:,ind(kk)),real_data_std(:,ind(kk)),'LineStyle','none','LineWidth',2)
    hold on
    
    TR = time_dep_TR_new_wSigma(nbd,final_energyi,TF_conc_t,RNAp_conc_t,mut_mat(ind(kk),:),vmax_array,group_array);
    plot(TF_conc_t,TR,'LineWidth',2)
    xlabel('TF concentration')
    ylabel('transcription rate')
    
    ylim([0 1000])
    title(string(title_name(ind(kk))))
end
saveas(gcf,append(file(1:end-4),'.jpeg'))
end
