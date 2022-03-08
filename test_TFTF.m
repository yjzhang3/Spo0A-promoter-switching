%% basic setup
TF_conc_t = @(x) x^3+30*x;
RNAp_conc = 5003;
nbd = 4;
tspan = 10;
Grnap = 10;

%% no TF-RNAp interaction (keep at 0), change TF-TF interaction
original_energy = [-1,7,2,Grnap,10,-20,0,-90,0,0];
pos_all = [5,6,8]; % same binding energy
% pos_all = [1,2,3,5,6,8]; % different binding energy

step = 1;
stepsize = 5;

[E_all,TR_all] = tuning(step,stepsize,pos_all,RNAp_conc,tspan,TF_conc_t,nbd,original_energy);

%% TF-RNAp interaction kept constant (but nonzero), change TF-TF interaction
original_energy = [-10,7,2,Grnap,10,-20,20,-90,-8,-10.4];
pos_all = [5,6,8]; % same binding energy
% pos_all = [1,2,3,5,6,8]; % different binding energy
step = 2;
stepsize = 2;

[E_all,TR_all] = tuning(step,stepsize,pos_all,RNAp_conc,tspan,TF_conc_t,nbd,original_energy);
%% no TF-TF interaction (keep at 0), change TF-RNAp interaction. Same binding energies
original_energy = [-10,3,2,Grnap,0,0,3,0,-8,-1];
pos_all = [7,9,10];

step = 2;
stepsize = 2;

[E_all,TR_all] = tuning(step,stepsize,pos_all,RNAp_conc,tspan,TF_conc_t,nbd,original_energy);

%% constant TF-TF interaction (nonzero), change TF-RNAp interaction. Same binding energies
% original_energy = [-10,78,2,Grnap,4,5,-1,6,-1,-1];
original_energy = [-10,300,2,Grnap,4,5,3,6,20,-1];
pos_all = [7,9,10];

step = 2;
stepsize = 2;

[E_all,TR_all] = tuning(step,stepsize,pos_all,RNAp_conc,tspan,TF_conc_t,nbd,original_energy);