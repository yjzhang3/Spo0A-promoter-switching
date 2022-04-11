%% a note on raw data
% Ps: column 2 and onward are one promoter Ps with all mutations. First
% column is true wild type
% Pv: column 2-9 are one promoter Pv with all mutations. First column is
% true wild type. Last column is when both promoters are deleted.
% all data start from T2 to T10 (index 1 to 9)
load('promoter_activity_single.mat')
load('Spo0A_dynamics.mat')
%% parameters
real_data_Ps = Ps_promoter_activity_mean(2:7,1:8); % exclude the real WT, instead WT is 

nbd = 4;
TF_conc_t = wtaCopy(1:6:end);
RNAp_conc = 0.0005;
mut_mat = [[1,1,1];[0,1,1];[1,0,1];[1,1,0];[0,0,1];[0,1,0];[0,0,1];[0,0,0]];
lb = zeros(10,1)-20;
ub = zeros(10,1)+20;
lb([5:6 8]) = -3;
ub([5:6 8]) = 3;
%%
pars = fit_data(nbd,TF_conc_t,RNAp_conc,mut_mat,real_data_Ps,lb,ub);

%% figures

% normalize the real data
for ll = 1:length(real_data_Ps(1,:))
    real_data_Ps(:,ll) = real_data_Ps(:,ll)./max(real_data_Ps(:,ll));
end

figure();
for kk = 1:length(mut_mat)
    subplot(4,2,kk)
    plot(TF_conc_t,real_data_Ps(:,kk),'LineWidth',4)
    ylim([0 1])
    xlabel('TF concentration')
    ylabel('transcription rate')
    title(sprintf('mutation %d',kk))
end

figure();
for yy = 1:length(mut_mat)
    subplot(4,2,yy)
    TR = time_dep_TR_new(nbd,pars(2,:),TF_conc_t,RNAp_conc,mut_mat(yy,:));
    plot(TF_conc_t,TR,'LineWidth',4)
    ylim([0 1])
    xlabel('TF concentration')
    ylabel('transcription rate')
    title(sprintf('mutation %d',yy))
end

%% now explore the probability of each configuration
config_all = dec2bin(0:(2^nbd-1), nbd) - '0';
tspan = length(TF_conc_t);
p_all_config = zeros(2^nbd,tspan);
energyi = pars(1,:);
mut = [1,1,1];

figure();
for cc = 1:length(config_all)
    subplot(4,4,cc)
    p_all = zeros(tspan,1);
    for rr = 1:tspan
        p_all(rr) = prob_per_config_new(nbd, config_all(cc,:),energyi,mut,TF_conc_t(rr),RNAp_conc);
    end
    
    p_all_config(cc,:) = p_all;
    plot(p_all,'LineWidth',4)
    ylim([0 1])
    xlabel('time')
    ylabel('probability')
    spec = sprintf('configuration %s',dec2bin(cc-1,4));
    title(spec)
    set(gca,'FontSize',11)
    [M,I] = max(p_all);
    spec1 = sprintf('time = %d',I);
    xline(I,'--',spec1)
    spec2 = sprintf('prob = %2f',M);
    yline(M,'--',spec2)

end
