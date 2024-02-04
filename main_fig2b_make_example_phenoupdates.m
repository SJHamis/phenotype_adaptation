%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% This code file produces example phenotype paths.
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
proffit_matrix=readmatrix('generate_PFM_fig1/PFM.xls');
dose_vec=readmatrix('generate_PFM_fig1/dosevector.xls');

% Chose drug
drug_c_in=300;
no_daily_cell_updates=2;

x_vec = linspace(0,1,11);
no_cells = 6;
no_data_cols =6;
end_t = 28;
x0_opt=0;

no_cells_init=no_cells;
cell_count_mat=zeros(6,end_t+1);
cell_count_mat(:,1)=no_cells_init*ones(6,1);
death_count_mat=zeros(6,end_t+1);
sp_idx=1;

% Define colour schemes.
mycol_candypop =[0 246 255
  170 0 255
  255 0 170
  255 170 0
  255 241 0
  170 255 0]/255;

tl=tiledlayout(3,4);
ft = 'Times';
fsz = 24; 

%Loop through 3 tratmrent options and 4 update options. 
%to: 1:no, 2:intermittent, 3:continous treatment.
%uo: 1:no, 2:diffusion, 3:advection-diffusion, 4:advection updates.
for to=1:3
    for uo=1:4
        treatment_opt=to;
        update_opt=uo;         
        no_cells = no_cells_init;
        no_data_cols =6;

        % Set the initial list of cells:
        cell_list = setFixedInitialCellList(no_cells,no_data_cols,proffit_matrix,x0_opt,dose_vec);
        time_evo = {cell_list};
        cell_counter = no_cells;


    
        % Progress in time:
        for t = 1:end_t+1

            %get drug concentration:
            drug_c = get_drug_c(t,treatment_opt,drug_c_in);

            %update phenotypic states:
            for k = 1:no_daily_cell_updates
                cell_list=updateCellFitness(cell_list,update_opt,drug_c,proffit_matrix,x0_opt,dose_vec);
            end
            
            new_cell_list=cell_list();
            death_counter=0;
    
            for c = 1:size(cell_list,1)
                if(cell_list(c,5)>0) %if positive growth rate
                    if(cell_list(c,5)>rand) %if division occurs 
                        new_cell_list = [new_cell_list; zeros(1,no_data_cols)]; %add daughter cell to list
                        cell_counter = cell_counter + 1;
                        new_cell_list(end,1) = cell_counter; %index
                        new_cell_list(end,2) = cell_list(c,2); %ancestor
                        new_cell_list(end,3:5) = cell_list(c,3:5); %set to be parent's fitness                    
                    end
                else %if negative growth rate
                    if(abs(cell_list(c,5))>rand)
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

        end

        % Plot the dynamics of phenotypic state x of all cells
        nexttile
        ofs=0.01; %offset used in plots
        % Make scheduling-shading
        if(to==2)
            fill([0 0 1+5*ofs 1+5*ofs],[0 28 28 0],'k','FaceAlpha',0.4,'LineStyle','none')
        elseif(to==3)
            fill([0 0 1+5*ofs 1+5*ofs],[0 7 7 0],'k','FaceAlpha',0.4,'LineStyle','none')
            hold on
            fill([0 0 1+5*ofs 1+5*ofs],[14 21 21 14],'k','FaceAlpha',0.4,'LineStyle','none')
        end
        hold on;

        max_cell_inx = time_evo{end}(end,1);    
        for idx = 1:max_cell_inx 
            t_vec=[];
            x_vec=[];

            %set empty ancestor
            anc = 0;
        
            for t=1:end_t+1
                if(ismember(idx,time_evo{t}(:,1)))
                    t_vec=[t_vec t];
                    tmp_row=find( time_evo{t}(:,1)==idx );
                    x_vec=[x_vec time_evo{t}(tmp_row,4)];

                    if anc ==0
                        anc = time_evo{1,t}(tmp_row,2);
                    end
                end
            end
            hold on;

            if(isempty(t_vec)==false)

                xoffset=(anc-1)*ofs;
                x_vec=x_vec+xoffset;
                t_vec=t_vec-1;

                pt=plot(x_vec,t_vec,'color',mycol_candypop(anc,:),'LineWidth',2);

                if(uo==1)
                %ylabel("days")
                end
                %xlabel("phenotype state (x)")
                set(gca, 'YDir','reverse','xaxisLocation','top')
                ylim([0 28])
                xlim([0 1+5*ofs])
                if(to==1)
                    xticks(linspace(0,1,no_cells_init));
                else
                    xticks([]);
                end

                yticks([0 7 14 21 28]);
                yticklabels({'0' '7', '14', '21', '28'});
        
            end
        end

        %%% prep data to plot cell count & death count dynamics
        cell_count_vec=zeros(1,end_t);
        for t=1:end_t
            cell_count_vec(t)=size(time_evo{t},1);
        end
        cell_count_mat(sp_idx,2:end)=cell_count_vec;

        set(findall(gcf,'type','text'), 'FontSize', fsz, 'Color', 'k','FontName', ft)
        set(gca,'FontSize', fsz, 'FontName', ft)

        sp_idx=sp_idx+1;
        %%%
    end
end

%fontsize(tl,16,"points")





