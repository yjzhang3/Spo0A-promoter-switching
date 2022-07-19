function TR = transcription_rate_new(nbd,energyi,mut,TF_conc,RNAp_conc,vmax)
% this is another way to incorporate mutations in the data
% since vmax is unique to strain only, given a particular mut, there should
% only be one vmax

%% exclude configurations that have neither promoters bound
bins = dec2bin(0:(2^nbd-1), nbd) - '0';
on_ind = find(bins(:,end) == 1); 
on_config = bins(on_ind,:);

% consider mutated sites (exclude those configurations where mutated site
% is 1)
final_config = on_config;
for ss = 1:length(mut) % ss range from 1 to 3
    if mut(ss) == 0
        final_ind = find(final_config(:,ss)~=1); % final_ind stores the index 
        % with respect to 8 configurations
        final_config = final_config(final_ind,:);
    end
end

final_ind = find(ismember(on_config,final_config,'row'));

if mut == [1,1,1]
    final_ind = 1:8;
end


TR = 0;
for oo = 1:length(final_config(:,1))
    ind = final_ind(oo); % now get the actual index with respect to 8 configurations
    if isvector(final_config)
        config = final_config;
    else
        config = final_config(oo,:);
    end
    curr_p = prob_per_config_new(nbd,config,energyi,mut,TF_conc,RNAp_conc);
    curr_tr = vmax*curr_p;
    TR = TR + curr_tr;
end

end