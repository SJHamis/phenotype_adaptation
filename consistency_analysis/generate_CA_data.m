%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Generate data shown in Figure 2 for consistency analysis (CA).
% Shows how phenotype distributions change over time in response
% to no (orange), 
% continous (purple) and 
% intermittent (green) treatments
% and the no/unbiased/semi-biased/biased update strategies (left to right)
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

%User inputs
% Set selected drug dose
drug_c_in=300;
% Set number of simulations (100 in preprint, 10 here for faster
% computation time)
no_reps = 100*20; %100 in manuscript
% Set number of cells at the simulation stat
no_cells = 11; %choose a multiple of 11 to make an evenly distributed start
no_daily_pheno_ups=2; %set number of phenotype updates per day

% No user inputs under this line

% Import the chosen growth rate matrix. 
proffit_matrix=readmatrix('generate_fitnessmatrix_and_fig1/fitnessmatrix.xls');
dose_vec=readmatrix('generate_fitnessmatrix_and_fig1/dosevector.xls');
x_vec = linspace(0,1,11);
no_data_cols = 6;
burnin = 28; % days before the assigned "day 1"
end_t = burnin + 28;

no_cells_init=no_cells;
cell_count_mat=zeros(no_data_cols,end_t+1);
cell_count_mat(:,1)=no_cells_init*ones(no_data_cols,1);
death_count_mat=zeros(no_data_cols,end_t+1);

sp_idx=1;
no_treat_opts = 3;
no_ud_opts = 4; 
no_tps = 6; 
no_x_states = 11;

gr_x0=0;gr_x1=0;

phenotype_distribution_matrix = zeros(no_treat_opts,no_ud_opts,no_reps,no_tps,no_x_states);

for to=1:3
    for uo=1:4
        [to uo]
        for rep=1:no_reps
            no_cells = no_cells_init;

            %Set uniform phenotype distribution:
            cell_list = setInitialCellList_uniform(no_cells,no_data_cols,proffit_matrix,dose_vec);
            time_evo = {cell_list};
            cell_counter = no_cells;
            
            %Start time
            t=1;
            for xval = 1:no_x_states
                phenotype_distribution_matrix(to,uo,rep,t,xval)=sum(cell_list(:,3)==xval);
            end
    
            %Loop days 2 to end
            for t = 2:end_t+1

                % Set daily drug dose
                if t<= burnin
                    drug_c=0;
                else
                    drug_c = get_drug_c(t-burnin,to,drug_c_in);
                end
                
                % Update phenotype states
                for updates = 1:no_daily_pheno_ups
                    cell_list=updateCellFitness(cell_list,uo,drug_c,proffit_matrix,dose_vec,gr_x0,gr_x1); 
                end
                new_cell_list=cell_list();
                               
                death_counter=0;
    
                for c = 1:size(cell_list,1)
                    if(cell_list(c,5)>0) %if positive growth
                        if(cell_list(c,5)>rand) %if division occurs, add daugther cell
                            new_cell_list = [new_cell_list; zeros(1,no_data_cols)];
                            cell_counter = cell_counter + 1;
                            new_cell_list(end,1) = cell_counter;%index
                            new_cell_list(end,2) = cell_list(c,2);%grandparent
                            new_cell_list(end,3:5) = cell_list(c,3:5);%set to be parent's fitnes                    
                        end
                    else %if negative growth rate
                        if(abs(cell_list(c,5))>rand) %if death occurs, remove cell
                            new_cell_list(c,6)=1;
                            death_counter=death_counter+1;
                        end
                    end
                end
    
                death_count_mat((to-1)*3+uo,t+1) = death_counter;
                %remove all cells with death-marker
                idx = (new_cell_list(:,6) == 1);
                new_cell_list(idx,:) = [];
                time_evo = [time_evo {new_cell_list}];
                cell_list = new_cell_list;


                for xval = 1:no_x_states
                    phenotype_distribution_matrix(to,uo,rep,t,xval)=sum(cell_list(:,3)==xval);
                end
    
            end
        end
    end
end
