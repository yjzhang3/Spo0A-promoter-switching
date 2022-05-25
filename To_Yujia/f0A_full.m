
function out = f0A_full(spo0A_p,muts,promoter,args, RNAP) 

% promoter = 'v' or 's'
% A: RNAP with sigma A
% H: RNAP with sigma H

t = num2cell(muts);
[box1,box2,box3] = deal(t{:});
t = num2cell(args);
[   k1,k2,k3,k12,k13,k23,k123, ...
    k0a,k1a,k2a,k3a,k12a,k13a,k23a,k123a, ...
    k0h,k1h,k2h,k3h,k12h,k13h,k23h,k123h] = deal(t{:});
if ~box1 
    k1=0;k12=0;k13=0;k123=0;
end
if ~box2 
    k2=0;k12=0;k23=0;k123=0;
end
if ~box3 
    k3=0;k23=0;k13=0;k123=0;
end

A=RNAP;H=RNAP;


b1=k1*spo0A_p;
b2=k2*spo0A_p;
b3=k3*spo0A_p;

b12=k12*spo0A_p^2;
b23=k23*spo0A_p^2;
b13=k13*spo0A_p^2;
b123=k123*spo0A_p^3;

b0a=A*k0a;
b1a=b1*A*k1a;
b2a=b2*A*k2a;
b3a=b3*A*k3a;
b12a=b12*A*k12a;
b13a=b13*A*k13a;
b23a=b23*A*k23a;
b123a=b123*A*k123a;

b0h=H*k0h;
b1h=b1*H*k1h;
b2h=b2*H*k2h;
b3h=b3*H*k3h;
b12h=b12*H*k12h;
b13h=b13*H*k13h;
b23h=b23*H*k23h;
b123h=b123*H*k123h;

if promoter=='v'
  Z=1+b1+b2+b3+b12+b13+b23+b123 ...
    +b0a+b1a+b2a+b3a+b12a+b13a+b23a+b123a;

  out=(b0a+b1a+b2a+b3a+b12a+b13a+b23a+b123a)/Z;

else
    if promoter=='s'
     Z=1+b1+b2+b3+b12+b13+b23+b123 ...
     +b0h+b1h+b2h+b3h+b12h+b13h+b23h+b123h;
    out=(b0h+b1h+b2h+b3h+b12h+b13h+b23h+b123h)/Z;
    end

end

