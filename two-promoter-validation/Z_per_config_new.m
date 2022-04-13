function Z = Z_per_config_new(config,energyi,TF_conc,type)

% energyi = [A,B,C,Pv,AB,AC,APv,BC,BPv,CPv,Ps,AB,AC,APs,BC,BPs,CPs]

% config = [A,B,C,Pv,Ps] 0 == unbound, 1 == bound
nbd = length(config); % only length(energyi)-1 sites are bound at one time 
% length(config)-2 is the number of TF motif sites

if type == 's' % Ps\
    energyi_final = energyi([1:length(config)-1 (nbd+nbd*(nbd-1)/2+1):end]);
end
if type == 'v'  % Pv or neither is  bound
    % (if no promoter is bound, it doesn't matter which energy array we
    % use  because the standard energy and interactive energy for anything
    % that involves promoter will be 0.
    energyi_final = energyi(1:(nbd+nbd*(nbd-1)/2));
end

G0 = stand_energy(config,energyi_final);

[num_arr,conc_arr] = gen_config_prof(config,TF_conc);

num_TF = num_arr(1);
num_RNAp = num_arr(end);

conc_TF = conc_arr(1);
conc_RNAp = conc_arr(end);

Z = exp(-G0)*(conc_TF^num_TF)*(conc_RNAp^num_RNAp);

end