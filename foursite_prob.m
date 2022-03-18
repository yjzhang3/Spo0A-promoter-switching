energyi = [10,10,10,10,0,0,-8,0,8,0];
mut = [1,1,1];
nbd = 4;
RNAp_conc = 500;
TF_conc_t = @(x) x^3+30*x;

%% first replicate the transcription rate graph

for tt = 1:500 
    p_on = transcription_rate_new(nbd,energyi,mut,TF_conc_t(tt),RNAp_conc);
    p_all(tt) = p_on;
end 

plot(p_all,'LineWidth',5)
xlabel('time')
ylabel('transcription rate')
set(gca,'FontSize',12)

%% now explore the probability of each configuration
config_all = [[0,0,0,1];[0,0,1,1];[0,1,0,1];[0,1,1,1];[1,0,0,1];[1,0,1,1];[1,1,0,1];[1,1,1,1]];

figure();
for cc = 1:length(config_all)
    subplot(2,4,cc)
    for rr = 1:500
        p_all(rr) = prob_per_config_new(nbd, config_all(cc,:), energyi,mut,TF_conc_t(rr),RNAp_conc);
    end
    plot(p_all,'LineWidth',4)
    xlabel('time')
    ylabel('probability')
    spec = sprintf('configuration %d',cc);
    title(spec)
    set(gca,'FontSize',11)
end 
    