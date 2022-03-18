function p_on = transcription_rate_new(nbd,energyi,mut,TF_conc,RNAp_conc)

bins = dec2bin(0:(2^nbd-1), nbd) - '0';
on_config = bins(bins(:,end)==1,:);

p_on = 0;
for oo = 1:length(on_config)
    p_on = p_on + prob_per_config_new(nbd,on_config(oo,:),energyi,mut,TF_conc,RNAp_conc);
end