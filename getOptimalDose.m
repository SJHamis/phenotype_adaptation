function opt_dose_idx = getOptimalDose(phenotypic_list,proffit_matrix,max_dose_idx)
    % zero drug sum
    sum = 0; 
    opt_dose_idx = 1;
    for c = 1:length(phenotypic_list)
        sum = sum  + proffit_matrix(1,phenotypic_list(c));
    end

    %9: up to 300
    %10: up to 1000
    for d = 2:max_dose_idx
        sum_new = 0;
        for c = 1:length(phenotypic_list)
            sum_new = sum_new + proffit_matrix(d,phenotypic_list(c));
        end
        if(sum_new<sum)
            sum = sum_new;
            opt_dose_idx = d;
        end
    end
end


