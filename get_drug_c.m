function drug_c = get_drug_c(t,treatment_opt,drug_c_in)

    if(treatment_opt==1) %No treatment
        drug_c=0;
    elseif(treatment_opt==2) %Continous treatment
        drug_c=drug_c_in;
    elseif(treatment_opt==3) %Intermittent (1 week on-1 week off).
        if(t<8)
            drug_c=drug_c_in;
        elseif(t<15)
            drug_c=0;
        elseif(t<22)
            drug_c=drug_c_in;  
        else
            drug_c=0;
        end
    end  

end
