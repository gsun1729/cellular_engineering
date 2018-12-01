function quotient = rxn_number(rm)
    r1 = rand;
    sum_R = 0;
    quotient = 0;
    for reaction_quotient = 1:length(rm)
        sum_R = sum_R + rm(reaction_quotient);
        if r1 < sum_R
            quotient = reaction_quotient;
            break
        end
    end
    return
end