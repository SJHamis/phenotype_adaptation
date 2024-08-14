close all
clear all

% Prepare figure for svg saving.
fig = figure; 
fig.Renderer = 'painter';
fig.Position = [0 0 800 200];
tl=tiledlayout(2,2);
tl.TileSpacing = 'loose';
ft = 'Arial';
fsz = 16; 

%load('workspace_500nM_costs_to_5_6')
fitness_matrix=readmatrix('generate_fitnessmatrix_and_fig1/fitnessmatrix.xls');
dose_vec=readmatrix('generate_fitnessmatrix_and_fig1/dosevector.xls');

lime_scheme=[240,255,0
201,255,0
107 142 35%163,255,0 
86,255,0
35,101,51]/255;

no_reps=100;
end_t=32;
sp=1;
no_cells_init=100;
no_data_cols=6;
cell_count_mat_int=zeros(no_reps,end_t);
ic_dir=0;

for T_lenght = [1 4]
    subplot(1,2,sp)
    sp=sp+1;
    for speed = [1 4]
        for uo=4    
            hold on
            u_on=speed;
            u_off=speed;
            gr_x0=-0.2866; %interpolated
            gr_x1=0.0741;
            drug_c=0;

            for to=3 %intermitt 
                for rep=1:no_reps
                    %initialise cell population
                    no_cells = no_cells_init;
                    cell_list = setInitialCellList_custom(no_cells,no_data_cols,fitness_matrix,dose_vec,uo,gr_x0,gr_x1,ic_dir);
                    cell_counter = no_cells; 
                    deathcount_vec=zeros(1,end_t);

                    %loop through time points
                    for t = 1:end_t
                        
                        %set drug        
                        if(mod(t-1,T_lenght)==0)
                            if(drug_c==500)
                                drug_c=0;
                            else
                                drug_c=500;
                            end
                        end

                        %Update cells fitness
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
    
                cell_count_mat_int(rep,1:end)=cell_count_vec;  %This increases in size
                end %end of reps                      
            end%end of to

        %mean of normalised counts, then count cose
        norm_int = mean(cell_count_mat_int./cell_count_mat_int(:,1));
        std_int = std(cell_count_mat_int./cell_count_mat_int(:,1));
    
        tv = 1:end_t;
        tv = [tv, fliplr(tv)];
        hold on
        fill(tv,[norm_int+std_int, fliplr(norm_int-std_int)],lime_scheme(3,:),...
            'LineStyle','none','FaceAlpha',0.2)
        hold on

        if(speed==1)

            no_rows = end_t/T_lenght;
            intervals = 1:T_lenght:end_t;
            intervals = reshape(intervals, [2,no_rows/2])';
            intervals = intervals-1;

            % Draw the grey boxes
            for k = 1:size(intervals, 1)
                x_box = [intervals(k, 1) intervals(k, 1) intervals(k, 2) intervals(k, 2)];
                y_box = [-1 1 1 -1] * 20; % Extend box height to cover the plot area
                patch(x_box, y_box, [0.9 0.9 0.9], 'EdgeColor', 'none', 'FaceAlpha', 0.3); % Grey color with transparency
            end

            hold on
            plot(1:end_t,mean(cell_count_mat_int./cell_count_mat_int(:,1)),':','color',lime_scheme(3,:),'LineWidth',3);
            else
                hold on
                plot(1:end_t,mean(cell_count_mat_int./cell_count_mat_int(:,1)),'color',lime_scheme(3,:),'LineWidth',3);
        end
        end

        xticks([1 7 14 21 28])
        xlim([1 28])
        ylim([0 5])
        ylabel('normalised cell count','FontSize',16)
        xlabel('time (days)','FontSize',16)

    end
end

saveas(fig, 'Figure3c.svg'); 

