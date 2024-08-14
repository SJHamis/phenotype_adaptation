%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% 
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all
close all

% initialise cell population
% column 1 = cell number (by instansiation)
% column 2 = parent cell (0 if no parent)
% column 3 = x idx
% column 4 = x value (the drug resistance)
% column 5 = pf proliferative fitness
% column 6 = death marker

% Import the chosen PFM. 
fitness_matrix=readmatrix('generate_fitnessmatrix_and_fig1/fitnessmatrix.xls');
dose_vec=readmatrix('generate_fitnessmatrix_and_fig1/dosevector.xls');



%daily_pheno_ups_list = [0 -14 -7 -2 1:10]% 1:1:10]; %-N means every N:th day (written this way for alorightm algebra)

daily_pheno_ups_list = 1:10;
delta_treat_list =1:10;
max_delta_treat = length(delta_treat_list);

no_reps =10;
x_vec = linspace(0,1,11);
ft = 'Arial';
fsz = 10;   
end_t = 28*2;


g2=-0.2866; %log lin interpol
g4=0.0741;
gr_x0 = g2;
gr_x1 = g4;

for ic_dir = 0:1
figure
sp=1;


%make empty matrices to populate with results
end_counts = zeros(3,length(daily_pheno_ups_list),max_delta_treat); 
avg_counts = zeros(3,length(daily_pheno_ups_list),max_delta_treat); 
avg_rates = zeros(3,length(daily_pheno_ups_list),max_delta_treat); 
avg_rates_bin = zeros(3,length(daily_pheno_ups_list),max_delta_treat); 
elim_time = zeros(3,length(daily_pheno_ups_list),max_delta_treat); 


for drug_dose_on = [300 500]
for update_opt=4:4
for upd_idx = 1:length(daily_pheno_ups_list)
    upd_idx
    for treat_idx = 1:length(delta_treat_list)
        for sim = 1:no_reps  

            delta_treat=delta_treat_list(treat_idx);
            no_daily_pheno_ups = daily_pheno_ups_list(upd_idx);
            trackmsg = [no_daily_pheno_ups,delta_treat,sim];

            %initiate cell population
            
            no_cells =10;
            no_data_cols = 6;            
            x0_opt = 0;
            
            drug_c_in = 0;
            %cell_list = setInitialCellList(no_cells,no_data_cols,fitness_matrix,x0_opt,dose_vec);
            cell_list = setInitialCellList_custom(no_cells,no_data_cols,fitness_matrix,dose_vec,update_opt,gr_x0,gr_x1,ic_dir);
            time_evo = {cell_list};
            cell_counter = no_cells;   
            
            dose_used = zeros(1,end_t+1);
            dose_used(1)=0;
            mycol = jet(10);        
            drug_idx = 1;
            drug_c = 0;
            sum_rate_tmp = 0; 
            rate_counter = 0;  
            tmp_elim_time_fastest=end_t+1;

            % For each time
            for t = 1:end_t

                %The drug is swithced when mod(t,delta_treat==0)
                if(treat_idx==length(delta_treat_list)+1) %if 28day
                    drug_c=drug_dose_on;
                else
                    if(mod(t-1,delta_treat)==0) %switch drug on mod0 days
                        if(drug_c==0)
                            drug_c=drug_dose_on;
                        else
                            drug_c=0;
                        end
                    end
                end

                dose_used(t+1)=drug_c; %add to drug tracker

                order_of_ops = rand();

                if(order_of_ops>0.5)

                    %update cells
                    if(no_daily_pheno_ups>=1)
                        for updates = 1:no_daily_pheno_ups %update the cells as many times as should
                            cell_list=updateCellFitness(cell_list,update_opt,drug_c,fitness_matrix,dose_vec,gr_x0,gr_x1);
                        end
                    elseif(no_daily_pheno_ups==0)
                        cell_list=updateCellFitness(cell_list,1,drug_c,fitness_matrix,dose_vec,gr_x0,gr_x1); %1 indicates no x upd (only d)
                    elseif(no_daily_pheno_ups<0)
                        tmp_T = -no_daily_pheno_ups;
                        if(mod(t,tmp_T)==0)
                            cell_list=updateCellFitness(cell_list,update_opt,drug_c,fitness_matrix,dose_vec,gr_x0,gr_x1);
                        else
                            cell_list=updateCellFitness(cell_list,1,drug_c,fitness_matrix,dose_vec,gr_x0,gr_x1); %1 indicates no x upd (only d)
                        end
                    end
                    %resolve death birth events
                    new_cell_list=cell_list;
                    death_counter=0;


                    %record rates
                    if(isempty(cell_list))
                        if(t<tmp_elim_time_fastest)
                            tmp_elim_time_fastest = t;
                        end
                    else
                        sum_rate_tmp = sum_rate_tmp + mean(cell_list(:,5));
                        rate_counter = rate_counter+1;
                    end

                    mean(cell_list(:,5));

                 

                    %resolve cell division & death
                    for c = 1:size(cell_list,1)
                        if(cell_list(c,5)>0) %if positive growth
                            if(cell_list(c,5)>rand) %if division occurs, %add daugther cell
                                new_cell_list = [new_cell_list; zeros(1,no_data_cols)];
                                cell_counter = cell_counter + 1;
                                new_cell_list(end,1) = cell_counter;%index
                                %new_cell_list(end,2) = cell_list(c,1);%parent
                                new_cell_list(end,2) = cell_list(c,2);%grandparent
                                new_cell_list(end,3:5) = cell_list(c,3:5);%set to be parent's fitnes
                            end
                        else
                            if(abs(cell_list(c,5))>rand)
                                new_cell_list(c,6)=1;
                                death_counter=death_counter+1;
                            end
                        end
                    end
                    %remove all cells with death-marker
                    idx = (new_cell_list(:,6) == 1);
                    new_cell_list(idx,:) = [];
                    time_evo = [time_evo {new_cell_list}];
                    cell_list = new_cell_list;
                
                
                else

                    
                    %update cells fitness through drug but nod phenotype
                    cell_list=updateCellFitness(cell_list,1,drug_c,fitness_matrix,dose_vec,gr_x0,gr_x1);


                    %record rates
                    if(isempty(cell_list))
                        if(t<tmp_elim_time_fastest)
                            tmp_elim_time_fastest = t;
                        end
                    else
                        sum_rate_tmp = sum_rate_tmp + mean(cell_list(:,5)) ;
                        rate_counter = rate_counter+1;
                    end



                    mean(cell_list(:,5));


                    %resolve death birth events
                    new_cell_list=cell_list;
                    death_counter=0;
                    %resolve cell division & death
                    for c = 1:size(cell_list,1)
                        if(cell_list(c,5)>0) %if positive growth
                            if(cell_list(c,5)>rand) %if division occurs, %add daugther cell
                                new_cell_list = [new_cell_list; zeros(1,no_data_cols)];
                                cell_counter = cell_counter + 1;
                                new_cell_list(end,1) = cell_counter;%index
                                %new_cell_list(end,2) = cell_list(c,1);%parent
                                new_cell_list(end,2) = cell_list(c,2);%grandparent
                                new_cell_list(end,3:5) = cell_list(c,3:5);%set to be parent's fitnes
                            end
                        else
                            if(abs(cell_list(c,5))>rand)
                                new_cell_list(c,6)=1;
                                death_counter=death_counter+1;
                            end
                        end
                    end
                    %remove all cells with death-marker
                    idx = (new_cell_list(:,6) == 1);
                    new_cell_list(idx,:) = [];
                    time_evo = [time_evo {new_cell_list}];
                    cell_list = new_cell_list;


                    %update cells
                    if(no_daily_pheno_ups>=1)
                        for updates = 1:no_daily_pheno_ups %update the cells as many times as should
                            cell_list=updateCellFitness(cell_list,update_opt,drug_c,fitness_matrix,dose_vec,gr_x0,gr_x1);
                        end
                    elseif(no_daily_pheno_ups==0)
                        cell_list=updateCellFitness(cell_list,1,drug_c,fitness_matrix,dose_vec,gr_x0,gr_x1); %1 indicates no x upd (only d)
                    elseif(no_daily_pheno_ups<0)
                        tmp_T = -no_daily_pheno_ups;
                        if(mod(t,tmp_T)==0)
                            cell_list=updateCellFitness(cell_list,update_opt,drug_c,fitness_matrix,x0_opt,dose_vec,gr_x0,gr_x1);
                        else
                            cell_list=updateCellFitness(cell_list,1,drug_c,fitness_matrix,dose_vec,gr_x0,gr_x1); %1 indicates no x upd (only d)
                        end
                    end



                end


          


            
            end %end of time loop
    
            %v is cell counts
            v = [];
            for t = 1:end_t+1
                v=[v length(time_evo{t})];
            end

            end_counts(sp,upd_idx,treat_idx) = ...
                end_counts(sp,upd_idx,treat_idx) + v(end)/no_reps;

            avg_counts(sp,upd_idx,treat_idx) = ...
                avg_counts(sp,upd_idx,treat_idx) + mean(v)/no_reps;

            avg_rates(sp,upd_idx,treat_idx) = ...
                avg_rates(sp,upd_idx,treat_idx) + (sum_rate_tmp/rate_counter)/no_reps;


            elim_time(sp,upd_idx,treat_idx) = tmp_elim_time_fastest;

        end % end of repl. sim
    end %end of delta tret
end %end of no daily pheno ops


if(drug_dose_on==500)
   
    target_r = gr_x1;

    elseif(drug_dose_on==300)
    target_r = fitness_matrix(9,11);
    elseif(drug_dose_on==1000)
    target_r = fitness_matrix(10,11);
end

for upd_idx = 1:length(daily_pheno_ups_list)
    for treat_idx = 1:length(daily_pheno_ups_list)
        tmp_r = avg_rates(sp,upd_idx,treat_idx);
        if(tmp_r<target_r)
            avg_rates_bin(sp,upd_idx,treat_idx)=1;
        end
    end
end

subplot(2,3,sp)
tmp_map = avg_rates(sp,:,:);
tmp_map=reshape(tmp_map,length(daily_pheno_ups_list),max_delta_treat);
heatmap(tmp_map)

subplot(2,3,3+sp)
tmp_map = avg_rates_bin(sp,:,:);
tmp_map=reshape(tmp_map,length(daily_pheno_ups_list),max_delta_treat);
heatmap(tmp_map)

sp=sp+1;
    
colormap(hot)
title('set schedule')
end
   

end
tmp_title=strcat('avg_rates_56h_ic',num2str(ic_dir),'.mat');
save(tmp_title,'avg_rates');
end

