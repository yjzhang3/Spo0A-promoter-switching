function TR = transcription_rate_new(nbd,energyi,mut,TF_conc,RNAp_conc,vmax)

% vmax is now an array
% create a map that maps a unique vmax to each configuration
bins = dec2bin(0:(2^nbd-1), nbd) - '0';
on_config = bins(bins(:,end)==1,:);

for gg = 1:length(on_config)
    str = string(on_config(gg,:));
    curr_config = append(str(1),str(2),str(3));
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
    curr_config = convertStringsToChars(append(str(1),str(2),str(3)));
    vmax = M(curr_config);
    curr_tr = vmax*prob_per_config_new(nbd,on_config(oo,:),energyi,mut,TF_conc,RNAp_conc);
    TR = TR + curr_tr;
end

TR = TR/(2^(nbd-1-sum(mut)));

end