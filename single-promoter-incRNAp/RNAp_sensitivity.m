syms x
f = @(x) x/(x+1);
df = @(x) 1/(x+1)^2;
sens = @(x) 1/(x+1);

tspan = 1:9.9;
for tt = 1:length(tspan)
    Yt(tt) = sens(tt);
end


