function TR = transcription_rate_new(nbd,energyi,mut,TF_conc,RNApH_conc,RNApA_conc,vmax)

% vmax is now an array
% create a map that maps a unique vmax to each configuration
bins = dec2bin(0:(2^nbd-1), nbd) - '0';
test = sum(bins(:,4:5)');
on_ind = find(test > 0);
on_config = bins(on_ind,:); % now we don't know what are
% considered "ON"
% so we will add up all probabilities with their vmax terms


%% we will first use a constant vmax to see the overall probability of each configuration
% if any configuration is unlikely to exist, then give that configuration
% vmax = 0

for gg = 1:length(on_config)
    str = string(on_config(gg,:));
    curr_config = append(str(1),str(2),str(3),str(4),str(5));
    keyset{gg} = convertStringsToChars(curr_config);
end

M = containers.Map(keyset,vmax);

% consider mutated sites
for ss = 1:length(mut)
    if mut(ss) == 0
        on_config(:,ss) = 0;
    end
end

TR = 0;
for oo = 1:length(on_config)
    str = string(on_config(oo,:));
    curr_config = append(str(1),str(2),str(3),str(4),str(5));
    curr_conf = convertStringsToChars(curr_config);
    curr_vmax = M(curr_conf);
    curr_tr = curr_vmax*prob_per_config_new(nbd,on_config(oo,:),energyi,mut,TF_conc,RNApH_conc,RNApA_conc);
%     curr_tr = prob_per_config_new(nbd,on_config(oo,:),energyi,mut,TF_conc,RNApH_conc,RNApA_conc);
    TR = TR + curr_tr;
end

TR = TR/(2^(3-sum(mut)));
% TR = vmax*TR/(2^(3-sum(mut)));

end