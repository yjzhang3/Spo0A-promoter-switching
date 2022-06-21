%% overall transcription rate
clear all
nbd = 4;

% Pv best
% load('Jun8-Pvonly-2vmax-b1.mat')
% load('A_RNAP.mat','A_conc_t')
% mut = [1,1,1];
% RNAp_conc_t = A_conc_t;

% Pv best
% load('Jun7-Psonly-1vmax-b9.mat')
% load('H_RNAP.mat','H_conc_t')
% mut = [1,1,1];
% RNAp_conc_t = H_conc_t;

% energyi = final_energyi;

% random 1
% energyi = [8,0,8,0,0,5,-8,0,0,-8];
energyi = [8,0,8,0,0,-5,-8,0,0,-8];
mut = [1,0,1];
TF_conc_t = get_aps(2:10);
RNAp_conc_t = [0.1,1,1.5,2,3,4,7.7,9,15]*10;
vmax_array.g1 = 200;
group_array.g1 = [1:8];


time_dep_TR_new_wSigma(nbd,energyi,TF_conc_t,RNAp_conc_t,mut,vmax_array,group_array)
%% individual probabilities and transcription rate
bins = dec2bin(0:(2^nbd-1), nbd) - '0';
on_ind = find(bins(:,end) == 1); 
on_config = bins(on_ind,:);

final_config = on_config;
for ss = 1:length(mut) % ss range from 1 to 3
    if mut(ss) == 0
        final_ind = find(final_config(:,ss)~=1); % final_ind stores the index 
        % with respect to 8 configurations
        final_config = final_config(final_ind,:);
    end
end

final_ind = find(ismember(on_config,final_config,'row'));

if mut == [1,1,1]
    final_ind = 1:8;
end

all_TR = zeros(length(TF_conc_t),length(final_config(:,1)));

%%
n_conf = length(on_config(:,1)); % should be 8, because it's all the possible configurations with
% RNAP bound and 3 other sites combo

% vmax is now an array
% create a map that maps a unique vmax to each configuration

% lets make vmax assignment using a function instead
vmax_final = vmax_assign(n_conf,vmax_array,group_array);

%%
figure();
for ii = 1:length(final_config(:,1))
    ind = final_ind(ii); % now get the actual index with respect to 8 configurations
    curr_config = final_config(ii,:);
    subplot(4,1,ii)
    
    curr_p = zeros(length(TF_conc_t),1);
    
    for kk = 1:length(TF_conc_t)
        curr_p(kk) = prob_per_config_new(nbd, curr_config, energyi, mut,TF_conc_t(kk),RNAp_conc_t(kk));
    end
    
    curr_vmax = vmax_final(ind);
    
    curr_TR = curr_vmax*curr_p;
    all_TR(:,ii) = curr_TR;
    
    
    plot(TF_conc_t,curr_TR,'LineWidth',3)
    
    xlabel('[TF]')
    ylabel('TR')
    ylim([0 200])
    title(append(dec2bin(ind-1,3),'1'))
end