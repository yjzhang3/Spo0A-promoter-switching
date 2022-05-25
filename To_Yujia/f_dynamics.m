function out = f_dynamics(muts,args,promoter,Aps,RNAP)
% muts: mutations. [box1, box2, box3], where box= 0 or 1 representing the box
% exists or not

% args: parameters to fit

% promoter: "v" or "s"

% Aps: spo0A_p levels at times 2:10

ts=2:10;
t=num2cell(args);
grs=[0.8402    0.7697    0.6288    0.4207    0.2311    0.1231    0.0766    0.0607    0.0580];
out=[];


for t = ts
   po=f0A_full(Aps(t-1),muts,promoter,args,RNAP);
   out=[out po/(grs(t-1)+1)*fglobal(grs(t-1),0.4)];
end

end

