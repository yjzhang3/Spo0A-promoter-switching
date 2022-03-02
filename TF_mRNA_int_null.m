function Q = TF_mRNA_int_null(list_int,mode)
% input: a list of TF-mRNA interactions per *site*
% output: factor Q that accounts for total TF-mRNA interacton per *config*

    
if mode == 1 % multiplicative model, meaning each TF-mRNA interaction could
    % depend on each other (for example, they all have to interact with
    % the same subunit of the transcipritonal machinery)
    Q = prod(list_int, 'all' ); 
end
if mode == 2
    Q = sum(list_int); % additive model, meaaning each TF-mRNA interaction 
    % is independent of each other
end
if Q == 0  
    Q = 1;
end

end