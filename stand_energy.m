function G0 = stand_energy(config,energyi)
% input:
% listi: list of binding box binding energies (same order as config
% profile), including that for mRNA

G0 = prod(config.*energyi, 'all' );
% if listi(i) = 0, meaning nothing is bound, no energy contribution from
% there (times 0 is 0), thus we need element multiplication

end


