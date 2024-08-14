%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Generate Figure 2c.
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
no_reps = 10; %100 in manuscript
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

%Make pretty color schemes
tangerine_scheme=[240,128,0 
244,187,68 
236,88,0 
255,165,0 
139, 64, 0]/255; 

lavender_scheme=[191, 148, 228
221, 160, 221
251, 174, 210
238, 130, 238
115, 79, 150
]/255;

lime_scheme=[240,255,0
201,255,0
163,255,0
86,255,0
35,101,51]/255;

%Time points to plot distributions
tp0 = 1;
tp=[0 7 14 21 28]+burnin;
tp=[tp0 tp];

%Place holders for phenotype counts
bars_no=zeros(11,1); %no treatment 
bars_cont=zeros(11,1); %continous treatment
bars_int=zeros(11,1); %intermittent treatment
mean_hist=zeros(4,3,6,11); %4 update rules, 3 treatment opts, 6 times, 11 phenotype options

% Plot variables
spix=1;
pos_list=[];
daylist={'day -28','day 1', 'day 7', 'day 14', 'day 21', 'day 28'};
dl_idx=1;

% Prepare figure for svg saving.
fig = figure; 
fig.Renderer = 'painter';
fig.Position = [0 0 800 600];
tl=tiledlayout(2,2);
tl.TileSpacing = 'loose';
ft = 'Arial';
fsz = 16; 

for tindx = 1:6
    for uo=1:4
       for to=1:3
            subplot(length(tp),no_ud_opts*no_treat_opts,spix) 
            for rep=1:no_reps

                %calculate total cell count at time t
                tot_count=0;
                for d = 1:11
                    tot_count=tot_count+phenotype_distribution_matrix(to,uo,rep,tp(tindx),d);
                end
        
                for d = 1:11
                    bars_no(d)=phenotype_distribution_matrix(to,uo,rep,tp(tindx),d)/tot_count; %normalise dist. with resp. to total cell count
                    mean_hist(uo,to,tindx,d) = mean_hist(uo,to,tindx,d)+(bars_no(d))/no_reps; %add to mean counter                
                end
    
                colindx=mod(rep,5)+1;

                if(to==1)
                    colscheme = tangerine_scheme;
                elseif(to==2)
                    colscheme = lavender_scheme;
                else
                    colscheme = lime_scheme;
                end

                b1=bar(1:11,bars_no,'FaceColor',colscheme(colindx,:),'Facealpha',0.1,'edgecolor','none','BarWidth',1);
                hold on
                if(rep==no_reps)
                    v_mean=zeros(1,11);
                    for dd=1:11
                        v_mean(dd)=mean_hist(uo,to,tindx,dd);
                    end
                    b1m=bar(1:11,v_mean,'FaceColor','none','edgecolor',colscheme(5,:),'BarWidth',1,'LineWidth',2);                
                end
            end

            xlim([0.5 11.5])
    
            %The below code creates custom subplot layout
            pos_tmp = get(gca, 'Position');
            pos_list=[pos_list;pos_tmp];

            % If not top row (t=-28days)
            if(spix>no_treat_opts*no_ud_opts)
                   pos_list(end,2)=pos_list(spix-no_treat_opts*no_ud_opts,2)-pos_tmp(4);
            end

            % If 2nd column in an update-rule
            if(mod(spix,3)==2)
                   pos_list(end,1)=pos_list(spix-1,1)+pos_tmp(3);
            end

            % If 3rd column in an update-rule
            if(mod(spix,3)==0)
                   pos_list(end,1)=pos_list(spix-2,1)+2*pos_tmp(3);
            end

            % If bottom row (last time point)
            if(spix>no_treat_opts*no_ud_opts*(no_tps-1))
                xticks([1 11])
                xticklabels({'0' '1'})
            else
                xticks([])
            end

            % If left-most plot on a row
            if(mod(spix,no_ud_opts*no_treat_opts)==1)
                ylabel(daylist(dl_idx));
                dl_idx=dl_idx+1;
            end
                 
            yticks([])
            
            set(gca, 'Position',[pos_list(end,1) pos_list(end,2) pos_tmp(3) pos_tmp(4)],'FontSize',12)
            set(findall(gcf,'type','text'), 'FontSize', fsz, 'Color', 'k','FontName', ft)
            set(gca,'FontSize', fsz, 'FontName', ft)

            spix=spix+1; 
        end

    end
end

saveas(fig, 'Figure2c.svg'); 
