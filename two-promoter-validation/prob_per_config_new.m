function p = prob_per_config_new(nbd,config,energyi,mut,TF_conc,type)

Zall = 0;

% config = [A,B,C,Pv,Ps] 0 == unbound, 1 == bound
if type == 's' % if Ps is bound
    [~,Zall_s] = Z_all_config_new(nbd,energyi,mut,TF_conc);
    Zall = Zall_s;
end

if type == 'v' % if Pv is bound or if no promoter is bound
    [Zall_v,~] = Z_all_config_new(nbd,energyi,mut,TF_conc);
    Zall = Zall_v;
end
    
Z = Z_per_config_new(config,energyi,TF_conc,type);

p = Z/Zall;

end 

