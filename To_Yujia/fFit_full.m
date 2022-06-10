function [x,fval,exitflag,output]=fFit_full(p_data)

options = optimoptions('particleswarm','Display','iter','SwarmSize',50,'MaxStallIterations',10);
% [ k1,k2,k3,k12,k13,k23,k123, ...
%   k0a,k1a,k2a,k3a,k12a,k13a,k23a,k123a, ...
%   k0h,k1h,k2h,k3h,k12h,k13h,k23h,k123h,...
%   A,H,kscale]
lb=[0,0,0,0,0,0,0,...
    0,0,0,0,0,0,0,0,...
    0,0,0,0,0,0,0,0,...
0.1,0.1,10];

ub=[200,200,200,200,200,200,200,...
    200,200,200,200,200,200,200,200,...
    200,200,200,200,200,200,200,200,...
1,1,500];

problem.options = options;
problem.solver = 'particleswarm';
problem.objective = @(x)f_error(p_data,x);
problem.lb = lb;
problem.ub = ub;
problem.nvars=26;

[x,fval,exitflag,output]=particleswarm(problem);

end