function final = weighted_msd(fake,real)
diff = zeros(length(real),1);

parfor jj = 1:length(fake)
    d = (fake(jj)-real(jj))^2; %every *element* of fake and real data subtractino
    df = d/real(jj)^2;
    diff(jj) = df;
end

final = 1/length(fake)*sum(diff,'all');
end
