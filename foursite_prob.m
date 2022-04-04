% energyi = [-5,-5,-5,-5,0,0,20,0,-5,-5];
energyi = x;
mut = [1,1,1];
nbd = 4;
RNAp_conc = 50;
TF_conc_t = @(x)x^3+30*x;
tspan = 50;

%% first replicate the transcription rate graph

for tt = 1:tspan
    p_on = transcription_rate_new(nbd,energyi,mut,TF_conc_t(tt),RNAp_conc);
    p_on_1 = transcription_rate(nbd,energyi,TF_conc_t(tt),RNAp_conc,mut);
    p_all(tt) = p_on;
%     p_all_1(tt) = p_on_1;
end 

figure();
plot(p_all,'r','LineWidth',5)
% hold on
% plot(p_all_1,'g-o')
xlabel('time')
ylabel('transcription rate')
set(gca,'FontSize',12)





%% now explore the probability of each configuration
config_all = dec2bin(0:(2^nbd-1), nbd) - '0';
p_all_config = zeros(2^nbd,tspan);

figure();
for cc = 1:length(config_all)
    subplot(4,4,cc)
    p_all = zeros(tspan,1);
    for rr = 1:tspan
        p_all(rr) = prob_per_config_new(nbd, config_all(cc,:), energyi,mut,TF_conc_t(rr),RNAp_conc);
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

%%
% for rr = 1:length(p_all_config(1,:))
%     [~,I] = max(p_all_config(:,rr))
% end
% 
check = sum(p_all_config);

