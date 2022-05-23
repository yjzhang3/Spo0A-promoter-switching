function TR_1 = transcription_rate_new_1(nbd,energyi,mut,TF_conc,RNApH_conc,RNApA_conc,vmax)
% this is another way to incorporate mutations in the data

%% exclude configurations that have neither promoters bound
bins = dec2bin(0:(2^nbd-1), nbd) - '0';
test = sum(bins(:,4:5)');
on_ind = find(test > 0); 
on_config = bins(on_ind,:); 

% vmax is now an array
% create a map that maps a unique vmax to each configuration
for gg = 1:length(on_config)
    str = string(on_config(gg,:));
    curr_config = append(str(1),str(2),str(3),str(4),str(5));
    keyset{gg} = convertStringsToChars(curr_config);
end

M = containers.Map(keyset,vmax);

clear on_ind

% consider mutated sites (exclude those configurations where mutated site
% is 1)
for ss = 1:length(mut) % ss range from 1 to 3
    if mut(ss) == 0
        on_ind = find(on_config(:,ss)~=1);
        on_config = on_config(on_ind,:);
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

TR_1 = TR;

end