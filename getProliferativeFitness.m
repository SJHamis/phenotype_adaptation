function pf = getProliferativeFitness(drug_c,x_idx,fitspace,dose_vec,gr_x0,gr_x1)
   
    if(drug_c==500)
        %Linear interploation from raw data (see main_generate_PFM.m file)
        gr_500_x0_raw= gr_x0;
        gr_500_x1_raw= gr_x1;
        fit_vec = linspace(gr_500_x0_raw,gr_500_x1_raw,11);
        pf=fit_vec(x_idx);
    else
        if(drug_c==0)
           drug_c=0.01;
        end
        didx = find(dose_vec==drug_c);
        
        fit_vec = fitspace(didx,:);
        pf=fit_vec(x_idx);
    end
    
end