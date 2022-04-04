energyi = x(2:end);
mut = [1,1,1];
nbd = 4;
RNAp_conc = x(1);

%% first replicate the transcription rate graph

for tt = 1:length(TF_conc_t)
    p_on = transcription_rate_new(nbd,energyi,mut,TF_conc_t(tt),RNAp_conc);
    p_all(tt) = p_on;
end 

plot(p_all,'LineWidth',5)
xlabel('time')
ylabel('transcription rate')
set(gca,'FontSize',12)

%% now explore the probability of each configuration (when promoter is bound)
config_all = [[0,0,0,1];[0,0,1,1];[0,1,0,1];[0,1,1,1];[1,0,0,1];[1,0,1,1];[1,1,0,1];[1,1,1,1]];

figure();
for cc = 1:length(config_all)
    subplot(2,4,cc)
    for rr = 1:length(TF_conc_t)
        p_all(rr) = prob_per_config_new(nbd, config_all(cc,:), energyi,mut,TF_conc_t(rr),RNAp_conc);
    end
    plot(p_all,'LineWidth',4)
    ylim([0 1])
    xlabel('time')
    ylabel('probability')
    spec = sprintf('configuration %d',cc);
    title(spec)
    set(gca,'FontSize',11)
end 

%% now explore the probability of each configuration (when promoter is not bound)
config_all = [[0,0,0,1];[0,0,1,1];[0,1,0,1];[0,1,1,1];[1,0,0,1];[1,0,1,1];[1,1,0,1];[1,1,1,1]];
config_all(:,4) = 0;

figure();
for cc = 1:length(config_all)
    subplot(2,4,cc)
    for rr = 1:length(TF_conc_t)
        p_all(rr) = prob_per_config_new(nbd, config_all(cc,:), energyi,mut,TF_conc_t(rr),RNAp_conc);
    end
    plot(p_all,'LineWidth',4)
    ylim([0 1])
    xlabel('time')
    ylabel('probability')
    spec = sprintf('configuration %d',cc);
    title(spec)
    set(gca,'FontSize',11)
end 

