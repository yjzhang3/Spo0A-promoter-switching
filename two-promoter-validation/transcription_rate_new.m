function [p_on_v,p_on_s] = transcription_rate_new(nbd,energyi,mut,TF_conc)

% config = [A,B,C,Pv,Ps] 0 == unbound, 1 == bound

bins = dec2bin(0:(2^(nbd-1)-1), nbd-1) - '0'; % all cases when Ps is bound
bins(:,nbd) = zeros(2^(nbd-1),1);
bins(:,nbd+1) = zeros(2^(nbd-1),1)+1;

bins_2 = dec2bin(0:(2^(nbd-1)-1), nbd-1) - '0'; % all cases when Pv is bound
bins_2(:,nbd) = zeros(2^(nbd-1),1)+1;
bins_2(:,nbd+1) = zeros(2^(nbd-1),1);

bins_final = [bins;bins_2]; % all possible configurations 
%  (Pv  and  Ps  can't bind  simulatenously)

% % consider mutated sites
% for ss = 1:length(mut)
%     if mut(ss) == 0
%         bins(:,ss) = 0;
%     end
% end

bins_v = bins_final(:,1:nbd);
bins_s = bins_final(:,[1:nbd-1 end]);

on_config_v = bins_v(bins_v(:,4) == 1,:);
on_config_s = bins_s(bins_s(:,4) == 1,:);

p_on_v = 0;
p_on_s = 0;
for oo = 1:length(on_config_v)
    p_s = prob_per_config_new(nbd,on_config_s(oo,:),energyi,mut,TF_conc,'s');
    p_v = prob_per_config_new(nbd,on_config_v(oo,:),energyi,mut,TF_conc,'v');
    p_on_v = p_on_v + p_v;
    p_on_s = p_on_s + p_s;
end

% p_on = p_on/(2^(nbd-1-sum(mut)));

end