function p_on = transcription_rate_new(nbd,energyi,mut,TF_conc,RNAp_conc)

bins = dec2bin(0:(2^nbd-1), nbd) - '0';
on_config = bins(bins(:,end)==1,:);

% consider mutated sites
for ss = 1:length(mut)
    if mut(ss) == 0
        on_config(:,ss) = 0;
    end
end

p_on = 0;
for oo = 1:length(on_config)
    p = prob_per_config_new(nbd,on_config(oo,:),energyi,mut,TF_conc,RNAp_conc);
    p_on = p_on + p;
end

p_on = p_on/(2^(nbd-1)/2^(sum(mut)));

end