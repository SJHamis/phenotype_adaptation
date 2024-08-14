%%%%
% Run after main_fig2d_300nM.m
%%%%


load('workspace_500nM_costs_BRAFV600E_fine.mat')

% Optional remove (not included in manuscript)
% remove = cost_list(:,2)  < 2;  % 2nd col <2 (i.e. 0 -not done- or 1 -not enough-)
% cost_list(remove,:) = [] ;  % remove those rows

%plot Data
figure(1)   
for uo = 1:no_uo_ops
    subplot(2,4,uo+no_uo_ops)
    plot(data_time, mean_cont_data,'^','color',lavender_scheme(5,:),'LineWidth',2)
    hold on
    plot(data_time, mean_int_data,'v','color',lime_scheme(3,:),'LineWidth',2)
    hold on
    errorbar(data_time,mean_cont_data,SEM_cont,'color',lavender_scheme(5,:),"LineStyle","none",'LineWidth',2)
    hold on
    errorbar(data_time,mean_int_data,SEM_int,'color',lime_scheme(3,:),"LineStyle","none",'LineWidth',2)
    ylim([0 15])
end

for uo=1:no_uo_ops

     cost_list_tmp = cost_list(cost_list(:,1)==uo,:);
     [min_val,row] = min(cost_list_tmp(:,end));
     figure(1)

     hold on
     subplot(2,no_uo_ops,uo+no_uo_ops)

     u_on=cost_list_tmp(row,2);
     u_off=cost_list_tmp(row,3);
     gr_x0=cost_list_tmp(row,4);
     gr_x1=cost_list_tmp(row,5);

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
    std_cont = std(cell_count_mat_cont./cell_count_mat_cont(:,1));
    std_int = std(cell_count_mat_int./cell_count_mat_int(:,1));

    tv = 1:end_t;
    tv = [tv, fliplr(tv)];
    hold on
    fill(tv,[norm_cont+std_cont, fliplr(norm_cont-std_cont)],lavender_scheme(5,:),...
        'LineStyle','none','FaceAlpha',0.2)
    hold on
    fill(tv,[norm_int+std_int, fliplr(norm_int-std_int)],lime_scheme(3,:),...
        'LineStyle','none','FaceAlpha',0.2)
    hold on

    plot(1:end_t,mean(cell_count_mat_cont./cell_count_mat_cont(:,1)),'color',lavender_scheme(5,:),'LineWidth',2);
    hold on
    plot(1:end_t,mean(cell_count_mat_int./cell_count_mat_int(:,1)),'color',lime_scheme(3,:),'LineWidth',2);
    tmp_title = num2str([u_on u_off gr_x0 gr_x1]);

    xticks([1 7 14 21 28])
    xlim([1 28])
    xlabel('time (days)','FontSize',13)

    if(uo==1)
        ylabel('normalised count','FontSize',13)
    end
    yticks([])
    ylim([0 15])
    yyaxis right
    set( gca, 'YTick', [] )
    set( gca, 'YColor', 'k' )
    yticks([1 5 12])
    ylim([0 15])
    parameter_cobmination = num2str([u_on u_off gr_x0 gr_x1])


    %%%DEATH
    figure(2)
    subplot(2,no_uo_ops,uo+no_uo_ops)

    mean_cont = mean(death_count_mat_cont(1:no_reps,:) ./cell_count_mat_cont(1:no_reps,:));
    std_cont = std(death_count_mat_cont(1:no_reps,:) ./cell_count_mat_cont(1:no_reps,:));
    mean_int = mean(death_count_mat_int(1:no_reps,:) ./cell_count_mat_int(1:no_reps,:));
    std_int = std(death_count_mat_int(1:no_reps,:) ./cell_count_mat_int(1:no_reps,:));
    
    tv = 1:end_t;
    tv = [tv, fliplr(tv)];
    fill(tv,[mean_cont+std_cont, fliplr(mean_cont-std_cont)],lavender_scheme(5,:),...
        'LineStyle','none','FaceAlpha',0.2)
    hold on
    plot(1:end_t,mean_cont,'color',lavender_scheme(5,:),'LineWidth',2,'LineStyle',':');
    hold on

    fill(tv,[mean_int+std_int, fliplr(mean_int-std_int)],lime_scheme(3,:),...
        'LineStyle','none','FaceAlpha',0.2)
    hold on
    plot(1:end_t,mean_int,'color',lime_scheme(3,:),'LineWidth',2,'LineStyle',':');
    hold on 

    xticks([1 7 14 21 28])
    xlim([1 28])
    xlabel('time (days)','FontSize',13)

    if(uo==1)
        ylabel('deaths per cell','FontSize',13)
    end

    ylim([0 0.4])
    yticks([])
    yyaxis right
    set( gca, 'YTick', [] )
    set( gca, 'YColor', 'k' )
    ylim([0 0.4])
    yticks([0 0.2 0.4])
end

saveas(fig1, 'Figure2d.svg'); 
saveas(fig2, 'Figure2f.svg'); 


