energyi = [10,15,15,0,-9,9];
mut = [1,1];
nbd = 3;
RNAp_conc = 1000;
config_all = [[0,0,0];[0,1,0];[1,0,0];[1,1,0]];

p_all = zeros(30,1);
TF_conc_t = @(x) 20*x^2;

%%
figure();
subplot(2,2,1)
for rr = 1:30
    p_all(rr) = prob_per_config_new(nbd, config_all(1,:), energyi,mut,TF_conc_t(rr),RNAp_conc);
end
plot(p_all)
xlabel('time')
ylabel('probability of configuration 000')

subplot(2,2,2)
for rr = 1:30
    p_all(rr) = prob_per_config_new(nbd, config_all(2,:), energyi,mut,TF_conc_t(rr),RNAp_conc);
end
plot(p_all)
xlabel('time')
ylabel('probability of configuration 010')

subplot(2,2,3)
for rr = 1:30
    p_all(rr) = prob_per_config_new(nbd, config_all(3,:), energyi,mut,TF_conc_t(rr),RNAp_conc);
end
plot(p_all)
xlabel('time')
ylabel('probability of configuration 100')

subplot(2,2,4)
for rr = 1:30
    p_all(rr) = prob_per_config_new(nbd, config_all(4,:), energyi,mut,TF_conc_t(rr),RNAp_conc);
end
plot(p_all)
xlabel('time')
ylabel('probability of configuration 110')


