function [TR_all,B] = change_energy(change,pos,RNAp_conc,tspan,TF_conc_t,nbd,B)
% change: an array of new energy values at a specific position in the array
% pos: index of that position in the energy array (only one at a time!)
% original: original energy array

% set the specific position to different energy values
B(:,pos) = change;

for nn = 1:length(B(:,1))
    
    TR_curr = time_dep_TR(TF_conc_t,tspan,RNAp_conc,nbd,B(nn,:));
    TR_all(nn,:) = TR_curr;
    
end 

    


