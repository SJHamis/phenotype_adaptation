
%%%%%%%%%%%
%
% Fit model parameters to in vitro data for cell counts 
% in 500 nM continous and intermittent treatments
%
%%%%%%%%%%%

% initialise cell population
% column 1 = cell number (by instansiation)
% column 2 = parent cell (0 if no parent)
% column 3 = x idx
% column 4 = x value (the drug resistance)
% column 5 = pf proliferative fitness
% column 6 = death marker

close all
clear all

% Import data
fitness_matrix=readmatrix('generate_fitnessmatrix_and_fig1/fitnessmatrix.xls');
dose_vec=readmatrix('generate_fitnessmatrix_and_fig1/dosevector.xls');
cf = readmatrix('generate_fitnessmatrix_and_fig1/invitro_data/kavran_28day_treatments.csv');

% User inputs
%Set model constants
drug_vec=[500];
no_reps = 10; %100 in manuscript
end_t=28;
no_uo_ops = 4;
no_cells_init = 100;
x0_opt=0;

%Set model params to fit
no_updates_on = 1:10;%updates per day when the drug is on. 0:10 in manuscript.
no_updates_off = 1:10; %updates per day when the drug is off. 0:10 in manuscrip.
no_n = 10; %fineness of below parameters to test. 100 in manuscript, 10 here for computational speed
gr_500_x0 = linspace(fitness_matrix(9,1),fitness_matrix(10,1),no_n);
gr_500_x1 = linspace(fitness_matrix(9,end),fitness_matrix(10,end),no_n);

% Handle data
% Each row corresponds to cell counts for 1 replicate
data_time = [1,7:7:28];
data_cont=cf(14:19,4:end);
data_int=cf(20:25,4:end);

% Normalise data and set outliers to zero.
data_cont=data_cont./data_cont(:,1);
data_int=data_int./data_int(:,1);
data_cont=data_cont.*(1-isoutlier(data_cont));
%data_int=data_int.*(1-isoutlier(data_int)); %Kavran et al. did not remove
%these outliers, so we don't either

mean_cont_data = zeros(1,length(data_time));
mean_int_data = zeros(1,length(data_time));
SEM_cont = zeros(1,length(data_time));
SEM_int = zeros(1,length(data_time));

no_data_cols = 6;
ic_dir = 0;

for t = 1:length(data_time)
    tmp_vec=[];
    tmp_vec_int=[];
    for rep = 1:size(data_cont,1)
        if(data_cont(rep,t)>0)
            tmp_vec = [tmp_vec data_cont(rep,t)];
        end
        if(data_int(rep,t)>0)
            tmp_vec_int = [tmp_vec_int data_int(rep,t)];
        end
    end
    mean_cont_data(t)=mean(tmp_vec);
    SEM_cont(t)=std(tmp_vec)/sqrt(length(tmp_vec));
    mean_int_data(t)=mean(tmp_vec_int);
    SEM_int(t)=std(tmp_vec_int)/sqrt(length(tmp_vec_int));
end

data_row_idx=1;

cell_count_mat_cont=zeros(no_reps,end_t);
death_count_mat_cont=zeros(no_reps,end_t);

cell_count_mat_int=zeros(no_reps,end_t);
death_count_mat_int=zeros(no_reps,end_t);

%cols: update opt, u_on, u_off, gr_x0, gr_x1, cost
cost_list = zeros(no_uo_ops*length(no_updates_on)*length(no_updates_off)*length(gr_500_x0)*length(gr_500_x1),6);

c_idx = 1;

for u_on = no_updates_on
    u_on
    for u_off = no_updates_off
        for gr_x0 = gr_500_x0 
            for gr_x1 = gr_500_x1             
                 for uo=1:no_uo_ops  %update models PUP1-4
                     %current_parameter_combination = [uo u_on u_off gr_x0 gr_x1] %for seeing progress

                     if(uo==1 && (u_on>1 || u_off>1))
                         %do nothing (for uo 1 there is no pheno adaptation)
                     else

                     for to=2:3 %no, cont, int
                        for rep=1:no_reps

                            %initialise cell population
                            no_cells = no_cells_init;
                            cell_list = setInitialCellList_custom(no_cells,no_data_cols,fitness_matrix,dose_vec,uo,gr_x0,gr_x1,ic_dir);
                            cell_counter = no_cells; 
                            deathcount_vec=zeros(1,end_t);
            
                            %loop through time points
                            for t = 1:end_t
                                %set drug and update
                                drug_c = get_drug_c(t,to,drug_vec(1)); 

                                if(drug_c>0)
                                    for updates = 1:u_on %uppdate when drug is on  
                                        cell_list=updateCellFitness(cell_list,uo,drug_c,fitness_matrix,dose_vec,gr_x0,gr_x1);
                                    end
                                else
                                    for updates = 1:u_off %uppdate when drug is off  
                                        cell_list=updateCellFitness(cell_list,uo,drug_c,fitness_matrix,dose_vec,gr_x0,gr_x1);
                                    end
                                end
                                new_cell_list=cell_list();
                                death_counter=0;

                                %Resolve cell divison and death
                                for c = 1:size(cell_list,1)
                                    if(cell_list(c,5)>0) %if positive growth
                                        if(cell_list(c,5)>rand) %if division occurs, %add daugther cell
                                            new_cell_list = [new_cell_list; zeros(1,no_data_cols)];
                                            cell_counter = cell_counter + 1;
                                            new_cell_list(end,1) = cell_counter;%index
                                            new_cell_list(end,2) = cell_list(c,1);%parent
                                            new_cell_list(end,3:5) = cell_list(c,3:5);%set to be parent's fitnes                    
                                        end
                                    else
                                        if(abs(cell_list(c,5))>rand)
                                            new_cell_list(c,6)=1;
                                            death_counter=death_counter+1;
                                        end
                                    end
                                end
 
                                %count death events
                                deathcount_vec(t)=death_counter;
                                %remove all cells with death-marker
                                idx = (new_cell_list(:,6) == 1);
                                new_cell_list(idx,:) = [];
                                if(t==1)
                                    time_evo = {cell_list}; %at time 0
                                else
                                    time_evo = [time_evo {new_cell_list}];
                                end
                                cell_list = new_cell_list;
                            end%end time loop

                            %%% save count and death data
                            cell_count_vec=zeros(1,end_t);
                            for t=1:end_t
                                cell_count_vec(t)=size(time_evo{t},1); %populate cell count over time
                            end

                            if(to==2)
                                cell_count_mat_cont(rep,1:end)=cell_count_vec; 
                                death_count_mat_cont(rep,1:end)=deathcount_vec;
                            elseif(to==3)
                                cell_count_mat_int(rep,1:end)=cell_count_vec;  %This increases in size
                                death_count_mat_int(rep,1:end)=deathcount_vec;
                            end     
                        end %end of reps                      
                    end%end of to

                    %mean of normalised counts, then count cose
                    norm_cont = mean(cell_count_mat_cont./cell_count_mat_cont(:,1));
                    norm_int = mean(cell_count_mat_int./cell_count_mat_int(:,1));
                    

                    cost_list(c_idx,:) = [uo,u_on,u_off,gr_x0,gr_x1,...
                        rmse([norm_cont(data_time) norm_int(data_time)],[mean_cont_data mean_int_data]) ];
                    c_idx=c_idx+1;   

                     end %end else

                 end
            end
        end
    end
end

remove = cost_list(:,6)  == 0; 
cost_list(remove,:) = [] ;

save('workspace_500nM_costs_BRAFV600E_fine')
