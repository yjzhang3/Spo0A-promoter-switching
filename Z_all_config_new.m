function Zall = Z_all_config_new(nbd,energyi,mut,TF_conc,RNAp_conc)

bins = dec2bin(0:(2^nbd-1), nbd) - '0'; % all possible configurations

% consider mutated sites
for ss = 1:length(mut)
    if mut(ss) == 0
        bins(:,ss) = 0;
    end
end

Zall = 0;
for bb = 1:length(bins)
    config = bins(bb,:);
    Z = Z_per_config_new(config, energyi,TF_conc,RNAp_conc);
    Zall = Zall + Z;
end

end
