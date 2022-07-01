function [RNAp_conc_t,vmax_array,diff] = estimate_RNAp_no0A4(ii,ind,energyi,vmax_array,group_array,file)

% inputs:
% i = type of promoter
% ind = which strain to fit in order to get RNAP dynamics
% should also be able to freely detemrine the grouping of vmax and have
% known or unknown vmax value

%% this file is trying to find how RNAp (more specificlaly, sigma-H bound
% holonenzyme concentration increaeses with time), and then incorporate
% this inofmration back to other strain types

% Ps: column 2 and onward are one promoter Ps with all mutations. First
% column is true wild type (both promoters exist)
% Pv: column 2-9 are one promoter Pv with all mutations. First column is
% true wild type (both promoters intact). Last column is when both promoters are deleted.
% all data start from T2 to T10 (index 1 to 9)
% we provide spo0A dynamics from hour  3 to 8 (so index into 2 to 7)

% no 0A4 box for now
load('promoter_activity_single.mat')

real_data_Ps = Ps_promoter_activity_mean(:,2:9); % exclude the real WT, instead WT is Ps only, no mutations
real_data_Ps_std = Ps_promoter_activity_std(:,2:9);
real_data_Pv = Pv_promoter_activity_mean(:,2:9);
real_data_Pv_std = Pv_promoter_activity_std(:,2:9);

TF_conc_t = get_aps(2:10);

if ii == 'v'
    real_data = real_data_Pv;
    real_data_std = real_data_Pv_std;
    title_name = {'Pv','1*4*','2*4*','3*4*','124*','134*','234*','1234*'};
    titlen = 'Pv-sigmaA';
end
if ii == 's'
    real_data = real_data_Ps;
    real_data_std = real_data_Ps_std;
    title_name = {'Ps','1*','2*','3*','12*','13*','23*','123*'};
    titlen = 'Ps-sigmaH';
end

%% other parameters

% if ind == 8:
% modify this accrodign to the strain you are fitting to
% energyi_sigma = zeros(10,1); % energy arrays used for estimating RNAP concentration
% doesn't matter what these energies are (set promoter energy to be 0, so
% that the found [RNAp] parameter incorporate energy value as well

nbd = 4;
% select whatever data selected  to find time-dependent RNAP
% treat both energies and RNAp as unknowns

data = real_data(:,ind); % the strains used for optimization
mut_mat = [[1,1,1];[0,1,1];[1,0,1];[1,1,0];[0,0,1];[0,1,0];[1,0,0];[0,0,0]]; 
mut = mut_mat(ind,:);

%% find [RNAp] as functino of time first
[pars,diff] = fit_data_sigma(nbd,energyi,TF_conc_t,mut,data,vmax_array,group_array,file);


RNAp_conc_t = pars(1:length(TF_conc_t));
vmax_raw = pars(length(TF_conc_t)+1:end);

%%

% create the vmax_array here:
% I would always put known vmax as g1!!! So this can avoid confusion and
% don't worry about assigning order
fn = fieldnames(group_array);
fnv = fieldnames(vmax_array);

if isempty(fnv)
    ol = 0;
end
if ~isempty(fnv)
    ol = numel(fieldnames(vmax_array)); % original length of vmax array field
end


for k=1:numel(fn)
    curr_key = (fn{k});
    if ~isfield(vmax_array,curr_key)  % only if this group doesn't exist, then we have unknown vmax
        curr_value = vmax_raw(k-ol);
        vmax_array.(fn{k}) = curr_value; % assign unknown parameters to the corresponding vmax group
    end
    % of the vmax's to be optimized the same as the order of the group
end


%% plot to see the fitting results
figure();
for kk = 1:length(ind)
    subplot(4,2,kk)
    
    errorbar(TF_conc_t,real_data(:,ind(kk)),real_data_std(:,ind(kk)),'LineStyle','none','LineWidth',2)
    hold on
    
    TR = time_dep_TR_new_wSigma(nbd,energyi,TF_conc_t,RNAp_conc_t,mut_mat(ind(kk),:),vmax_array,group_array);
    plot(TF_conc_t,TR,'LineWidth',2)
    xlabel('TF concentration')
    ylabel('transcription rate')
    
    ylim([0 1000])
    title(string(title_name(ind(kk))))
end
saveas(gcf,append(file(1:end-4),'.jpeg'))
end
