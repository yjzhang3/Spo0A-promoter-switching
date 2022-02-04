function TR_overall = transcription_rate(tspan,nbd,tfunc,G0)
% input: 
% tspan: time for simulation
% nbd: # of binding sites
% tfunc: how TF increases over time
% G0 = standard free energy

% output: transcription rate

% generate time dependent configuration profiles
t_overall = gen_fake_tdata(tspan,nbd,tfunc);
conc_overall = zeros(tspan,1);

TR_overall = zeros(tspan,0);
for tt = 1:tspan
    % partition function Z_on
    lb = (tt-1)*2^nbd+1;
    ub = tt*2^nbd;
    Zon = Z_on(nbd,t_overall(lb:ub),G0);

    % partition function Z_off
    Zoff = Z_off(nbd,t_overall(lb:ub),G0);
    
    ca = t_overall(lb).conc_arr;
    conc_overall(tt) = ca(1);
    
    TR_overall(tt) = Zon/(Zon+Zoff);
end

plot(1:tspan,TR_overall)
hold on
plot(1:tspan,conc_overall)
xlabel('time')
legend('TR','conc of TF')