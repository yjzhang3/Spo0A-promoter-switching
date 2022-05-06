function TR = transcription_rate_new(nbd,energyi,mut,TF_conc,RNAp_conc,vmax)

% vmax is now an array

bins = dec2bin(0:(2^nbd-1), nbd) - '0';
on_config = bins(bins(:,end)==1,:);

% keyset = {[0,0,0,1] [0,0,1,1] [0,1,0,1] [0,1,1,1] [1,0,0,1] [1,0,1,1] ...
%     [1,1,0,1],[1,1,1,1]};


% consider mutated sites
for ss = 1:length(mut)
    if mut(ss) == 0
        on_config(:,ss) = 0;
    end
end

TR = 0;
for oo = 1:length(on_config)
    config_curr = on_config(oo,:);
    config_type = bi2de(config_curr(1:3)); % the decimal representation of the configuration
    curr = vmax(config_type+1)*prob_per_config_new(nbd,on_config(oo,:),energyi,mut,TF_conc,RNAp_conc);
    TR = TR + curr;
end

TR = TR/(2^(nbd-1-sum(mut)));

end