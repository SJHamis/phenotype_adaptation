% initialise cell population
% column 1 = cell number (by instansiation)
% column 2 = parent cell (0 if no parent)
% column 3 = x idx
% column 4 = x value (the drug resistance)
% column 5 = pf proliferative fitness
% column 6 = death marker

close all
clear all

fitness_matrix=readmatrix('generate_fitnessmatrix_and_fig1/fitnessmatrix.xls');
dose_vec=readmatrix('generate_fitnessmatrix_and_fig1/dosevector.xls');

no_daily_pheno_ups = 2;
no_cells_init = 100;%100 in manuscript
no_data_cols = 6;
end_t = 28;
no_uo_ops = 4;
no_reps=100;
ic_dir = 0; %cells startt in x=0 for the directed update strategies

drug_vec=[300];

output_data_cont = zeros(length(drug_vec)*no_uo_ops,end_t);
output_data_int = zeros(length(drug_vec)*no_uo_ops,end_t);
data_row_idx=1;

gr_x0=0;gr_x1=0;

%Prepare for svg figure save
fig1 = figure(1);  
fig1.Renderer = 'painter';
fig1.Position = [0 0 800 300];
tl=tiledlayout(2,2);
tl.TileSpacing = 'loose';
fig2 = figure(2);  
fig2.Renderer = 'painter';
fig2.Position = [0 0 800 300];
tl=tiledlayout(2,2);
tl.TileSpacing = 'loose';

ft = 'Arial';
fsz = 16; 

%Custom colors
gray_scheme = [238, 238, 238
        174,174,174
        144, 144, 144
        92,92,92
        56,56,56
        ]/255;

lavender_scheme=[191, 148, 228
        221, 160, 221
        251, 174, 210
        238, 130, 238
        115, 79, 150
        ]/255;

lime_scheme=[240,255,0
        201,255,0
        107 142 35%163,255,0 
        86,255,0
        0,255,30]/255;

for dose_sett = 1:length(drug_vec)
    for uo=1:no_uo_ops  
        
        drug_c_in=drug_vec(dose_sett);
        sp_idx=1; %for plotting
        cell_count_mat=zeros(length(drug_vec),end_t);
        death_count_mat=zeros(length(drug_vec),end_t);

        for to=2:3 %1:no, 2:cont, 3:int treatment
            for rep=1:no_reps

                treatment_opt=to;
                update_opt=uo; 
                no_cells = no_cells_init;
                cell_list = setInitialCellList_custom(no_cells,no_data_cols,fitness_matrix,dose_vec,uo,gr_x0,gr_x1,ic_dir);
                cell_counter = no_cells; 
                deathcount_vec=zeros(1,end_t);

                for t = 1:end_t
                    drug_c = get_drug_c(t,treatment_opt,drug_c_in);

                    for updates = 1:no_daily_pheno_ups
                            cell_list=updateCellFitness(cell_list,update_opt,drug_c,fitness_matrix,dose_vec,gr_x0,gr_x1);
                    end

                    new_cell_list=cell_list();
                    death_counter=0;

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

                end

                %%% save data to (later)plot cell count & death count dynamics
                cell_count_vec=zeros(1,end_t);
                for t=1:end_t
                    cell_count_vec(t)=size(time_evo{t},1);
                end
                cell_count_mat(sp_idx,1:end)=cell_count_vec;  
                death_count_mat(sp_idx,1:end)=deathcount_vec;

                sp_idx=sp_idx+1;
                %%%
            end
        end


        tvec=1:28;

        figure(1)
        set(gcf)

        if(dose_sett==1)
            subplot(length(drug_vec)+1,no_uo_ops,uo+(dose_sett-1)*no_uo_ops)
        else
            subplot(length(drug_vec)+1,no_uo_ops,uo+no_uo_ops*2)
        end

        mean_cont = mean(cell_count_mat(1:no_reps,:)./cell_count_mat(1:no_reps,1));
        std_cont = std(cell_count_mat(1:no_reps,:)./cell_count_mat(1:no_reps,1));
        mean_int = mean(cell_count_mat(no_reps+1:no_reps*2,:)./cell_count_mat(no_reps+1:no_reps*2,1));
        std_int = std(cell_count_mat(no_reps+1:no_reps*2,:)./cell_count_mat(no_reps+1:no_reps*2,1));
        tv = 1:end_t;
        tv = [tv, fliplr(tv)];
        fill(tv,[mean_cont+std_cont, fliplr(mean_cont-std_cont)],lavender_scheme(5,:),...
            'LineStyle','none','FaceAlpha',0.2)
        hold on
        plot(1:end_t,mean_cont,'color',lavender_scheme(5,:),'LineWidth',2);
        hold on

        fill(tv,[mean_int+std_int, fliplr(mean_int-std_int)],lime_scheme(3,:),...
            'LineStyle','none','FaceAlpha',0.2)
        hold on
        plot(1:end_t,mean_int,'color',lime_scheme(3,:),'LineWidth',2);
        hold on

        plot(1:end_t,5*ones(length(1:end_t)),':','color',gray_scheme(2,:),'LineWidth',1)

        xlim([1 28])
        xticks([1 7 14 21 28])

        if(uo==1)
            ylabel('normalised count','FontSize',13)
        end

        if(uo<=2)
            yticks([])
            ylim([0 7.5])
            yyaxis right
            set( gca, 'YTick', [] )
            set( gca, 'YColor', 'k' )
            yticks([1 5])
            ylim([0 7.5])
        else
            yticks([])
            ylim([0 50])
            yyaxis right
            set( gca, 'YTick', [] )
            set( gca, 'YColor', 'k' )
            yticks([1 5 25 50])
            ylim([0 50])
        end

            
        %Save data
        output_data_cont(data_row_idx,:) = mean_cont;
        output_data_int(data_row_idx,:) = mean_int;
        data_row_idx=data_row_idx+1;
        %%%%

        figure(2)
        set(gcf)
        if(dose_sett==1)
            subplot(length(drug_vec)+1,no_uo_ops,uo+(dose_sett-1)*no_uo_ops)
        else
            subplot(length(drug_vec)+1,no_uo_ops,uo+no_uo_ops*2)
        end

        cell_count_mat(cell_count_mat==0) = -1 ;

        mean_cont = mean(death_count_mat(1:no_reps,:) ./cell_count_mat(1:no_reps,:));
        std_cont = std(death_count_mat(1:no_reps,:) ./cell_count_mat(1:no_reps,:));
        mean_int = mean(death_count_mat(no_reps+1:no_reps*2,:) ./ cell_count_mat(no_reps+1:no_reps*2,:));
        std_int = std(death_count_mat(no_reps+1:no_reps*2,:) ./cell_count_mat(no_reps+1:no_reps*2,:));
              
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

        xticks([])
        xlim([1 28])
        xticks([1 7 14 21 28])
 
        if(uo==1)
            ylabel('deaths per cell','FontSize',13)
        end

        yticks([])
        ylim([0 0.2])
        yyaxis right
        set( gca, 'YTick', [] )
        set( gca, 'YColor', 'k' )
        yticks([0 0.1 0.2])
        ylim([0 0.2])

    end%end update opt
    hold on

    if(uo==4 && dose_sett==5)
    legend('','continous treatment','','intermittent treatment','Orientation','horizontal','Location','bestoutside')
    end

end%end dose opt

%Export as high resolution png to avoid error bars getting chopped from svg
%exportgraphics(fig1,'fig2d.png','Resolution',1000)
%exportgraphics(fig2,'fig2f.png','Resolution',1000)

