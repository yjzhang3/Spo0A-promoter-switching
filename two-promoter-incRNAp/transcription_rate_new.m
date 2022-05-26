function TR = transcription_rate_new(nbd,energyi,mut,TF_conc,RNApH_conc,RNApA_conc,vmax_array,group_array)
% this is another way to incorporate mutations in the data

%% exclude configurations that have neither promoters bound

bins = dec2bin(0:(2^nbd-1), nbd) - '0';
test = sum(bins(:,4:5)');
on_ind = find(test > 0);
on_config = bins(on_ind,:); % now we don't know what are
% considered "ON"
% so we will add up all probabilities with their vmax terms

% consider mutated sites (exclude those configurations where mutated site
% is 1)
final_config = on_config;
final_ind = [];
for ss = 1:length(mut) % ss range from 1 to 3
    if mut(ss) == 0
        final_ind = find(final_config(:,ss)~=1); % final_ind stores the index 
        % with respect to 8 configurations
        final_config = final_config(final_ind,:);
    end
end

if mut == [1,1,1]
    final_ind = 1:24;
end

n_conf = length(on_config(:,1)); % should be 24, because it's all the possible configurations with
% RNAP bound and 3 other sites combo

% vmax is now an array
% create a map that maps a unique vmax to each configuration

% lets make vmax assignment using a function instead
vmax_final = vmax_assign(n_conf,vmax_array,group_array);
% this should also have 8 terms, universal for all mutation types


TR = 0;
for oo = 1:length(final_config(:,1))
    ind = final_ind(oo); % now get the actual index with respect to 24 configurations
    config = final_config(oo,:);
    curr_vmax = vmax_final(ind);
    curr_tr = curr_vmax*prob_per_config_new(nbd,config,energyi,mut,TF_conc,RNApH_conc,RNApA_conc);
    TR = TR + curr_tr;
end

end