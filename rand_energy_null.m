function rande = rand_energy_null(original_energy,n,posk)
% nbd: number of binding sites (plus promoter)
% n: number of different energy arrays you want to test
% posc: position to keep (at these positions, energy stays the same)

rande = randi([-10, 1000],n,length(original_energy));
for ii = 1:length(posk)
    rande(:,posk(ii)) = original_energy(posk(ii));
    
end
