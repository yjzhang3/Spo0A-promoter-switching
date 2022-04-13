function [Zall_v,Zall_s] = Z_all_config_new(nbd,energyi,mut,TF_conc)

% energyi = [A,B,C,Pv,AB,AC,APv,BC,BPv,CPv,Ps,AB,AC,APs,BC,BPs,CPs]
% config = [A,B,C,Pv,Ps] 0 == unbound, 1 == bound
% nbd =  4 

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

Zall_v = 0;
for bb = 1:length(bins_v)
    config = bins_v(bb,:);
    Z = Z_per_config_new(config, energyi,TF_conc,'v');
    Zall_v = Zall_v + Z;
end

Zall_s = 0;
for bb = 1:length(bins_s)
    config = bins_s(bb,:);
    Z = Z_per_config_new(config, energyi,TF_conc,'s');
    Zall_s = Zall_s + Z;
end

% Zall = Zall/2^(nbd-1-sum(mut));
% divide by 2^(number of mutations), where number of mutations is 
% nbd-1-sum(mut), becauase nbd-1 is the number of binding sites besides the
% promoter, and sum(mut) is the number of non-mutated sites

% [1,1,1] all working
% [0,0,0] all mutated
end
