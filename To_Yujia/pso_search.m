% find n solutions
n=5;
full_sols={};
for i=1:n
[x,fval,exitflag,output]=fFit_full(p_data);
full_sols{i}={x,fval,exitflag,output};

end
save('full_sols.mat','full_sols');