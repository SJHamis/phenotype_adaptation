function new_cell_list= updateCellFitness(cell_list,update_opt,drug_c,fitspace,dose_vec,gr_x0,gr_x1)

    lambda = 1;
    nu = 1;

    for c = 1:size(cell_list,1)

            x_idx=cell_list(c,3);
            x_idx_new=x_idx;

            if(update_opt==1)
                 % do nothing

            elseif(update_opt==2)
                tmp_dice = rand;
                if(tmp_dice>(1-lambda/2) && x_idx<11)
                    x_idx_new=x_idx+1;
                else
                    if(x_idx>1)
                        x_idx_new=x_idx-1;
                    end
                end
    
            elseif(update_opt==3)
                x_idx_prop = 0;
                tmp_dice = rand;    
                if(tmp_dice>(1-lambda/2) && x_idx<11)
                    x_idx_prop=x_idx+1;
                else
                    if(x_idx>1)
                        x_idx_prop=x_idx-1;
                    end
                end
                if(x_idx_prop>0 && ...
                        getProliferativeFitness(drug_c,x_idx_prop,fitspace,dose_vec,gr_x0,gr_x1)...
                        >getProliferativeFitness(drug_c,x_idx,fitspace,dose_vec,gr_x0,gr_x1))
                   x_idx_new=x_idx_prop;
                end
    
            elseif(update_opt==4)
                tmp_dice = rand;
                if(nu>tmp_dice)
                    if(x_idx<11 && ...
                            getProliferativeFitness(drug_c,x_idx+1,fitspace,dose_vec,gr_x0,gr_x1)...
                            >getProliferativeFitness(drug_c,x_idx,fitspace,dose_vec,gr_x0,gr_x1))
                        x_idx_new = x_idx+1;
                    elseif(x_idx>1 && ...
                            getProliferativeFitness(drug_c,x_idx-1,fitspace,dose_vec,gr_x0,gr_x1)...
                            >getProliferativeFitness(drug_c,x_idx,fitspace,dose_vec,gr_x0,gr_x1))
                        x_idx_new = x_idx-1;
                    end
                end
        
            end

            cell_list(c,3)=x_idx_new;
            cell_list(c,4)=(cell_list(c,3)-1)/10;
            cell_list(c,5) = getProliferativeFitness(drug_c,cell_list(c,3),fitspace,dose_vec,gr_x0,gr_x1);
    end
    new_cell_list=cell_list();
end


