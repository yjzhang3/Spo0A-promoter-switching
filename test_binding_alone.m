TF_conc_t = @(x) x^3+30*x;
% RNAp_conc = 100000;
RNAp_conc = 2444;
nbd = 4;
tspan = 10;
% original = [-10,78,1000,10,0,0,0,0,0,0];
original = [-10,78,493,10,0,0,0,0,0,0];
pos_all = [1,2,3];

step = 1;
stepsize = 56;

[E_all,TR_all] = tuning(step,stepsize,pos_all,RNAp_conc,tspan,TF_conc_t,nbd,original);


% as expected, when there is no protein-protein interaction, transcription
% rate is sole dependent on the product of mRNA concentration and its
% binidng energy. TF binding energy does not matter regardless if they are
% equal or not