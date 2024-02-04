%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Make phenotype distibutions 
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

drug_c_in=300;
no_cells = 100;
no_daily_pheno_ups=2; %10 to fully push dif-adv. 2 to fully push adv.


x_vec = linspace(0,1,11);
no_data_cols =6;
end_t = 28;
x0_opt=0;
no_cells_init=no_cells;
cell_count_mat=zeros(6,end_t+1);
cell_count_mat(:,1)=no_cells_init*ones(6,1);
death_count_mat=zeros(6,end_t+1);


sp_idx=1;
no_treat_opts = 3;
no_ud_opts = 4;
no_reps = 10;
no_tps = 5;
no_x_states = 11;

phenotype_distribution_matrix = zeros(no_treat_opts,no_ud_opts,no_reps,no_tps,no_x_states);

for to=1:3
    for uo=1:4
        for rep=1:no_reps

            treatment_opt=to;
            update_opt=uo;
            no_cells = no_cells_init;
            no_data_cols =6;
            cell_list = setInitialCellList(no_cells,no_data_cols,proffit_matrix,x0_opt,dose_vec);
            time_evo = {cell_list};
            cell_counter = no_cells;

            t=1;
            for xval = 1:no_x_states
                phenotype_distribution_matrix(to,uo,rep,t,xval)=sum(cell_list(:,3)==xval);
            end
    
            for t = 2:end_t+1
                drug_c = get_drug_c(t,treatment_opt,drug_c_in);

                for updates = 1:no_daily_pheno_ups
                    cell_list=updateCellFitness(cell_list,update_opt,drug_c,proffit_matrix,x0_opt,dose_vec);             
                end
                new_cell_list=cell_list();
                               
                death_counter=0;
    
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
    
                death_count_mat((to-1)*3+uo,t+1) = death_counter;
                %remove all cells with death-marker
                idx = (new_cell_list(:,6) == 1);
                new_cell_list(idx,:) = [];
                time_evo = [time_evo {new_cell_list}];
                cell_list = new_cell_list;


                for xval = 1:no_x_states
                    phenotype_distribution_matrix(to,uo,rep,t,xval)=sum(cell_list(:,3)==xval);
                    %phenotype_distribution_matrix(to,uo,rep,t,xval)=sum(cell_list(:,3)==xval);
                end
    
            end
        end
    end
end

dumn=zeros(11,1); %no treat
dumc=zeros(11,1); %cont treat
dumi=zeros(11,1); %int treat

spix=1;

gray_scheme = [238, 238, 238
174,174,174
144, 144, 144
92,92,92
56,56,56
]/255;

tangerine_scheme=[240,128,0 %tan
244,187,68 %MAN
236,88,0 %PERS
255,165,0 %OR
255,117,24]/255; %PUMP

lavender_scheme=[115, 79, 150
191, 148, 228
221, 160, 221
251, 174, 210
238, 130, 238
]/255;

lime_scheme=[240,255,0
201,255,0
163,255,0
86,255,0
0,255,30]/255;


tp=[1 7 14 21 28];

pos_list=[];

tolist={ 'no', 'cont.', 'interm.'};
tl_idx=1;
daylist={'day 1', 'day 7', 'day 14', 'day 21', 'day 28'};
dl_idx=1;

figure
ft = 'Times';
fsz = 16; 

for tindx = 1:5
    for uo=1:4
       for to=1:3
            subplot(length(tp),no_ud_opts*no_treat_opts,spix)
            for rep=1:no_reps

                %phenotype_distribution_matrix(to,uo,rep,t,xval)
                for d = 1:11
                    dumn(d)=phenotype_distribution_matrix(1,uo,rep,tp(tindx),d);
                    dumc(d)=phenotype_distribution_matrix(2,uo,rep,tp(tindx),d);
                    dumi(d)=phenotype_distribution_matrix(3,uo,rep,tp(tindx),d);
                end
    
                colindx=mod(rep,5)+1;
    
                if(to==1)
                    b1=bar(1:11,dumn,'FaceColor',tangerine_scheme(colindx,:),'Facealpha',.1,'edgecolor','none','BarWidth',1);
                    hold on
                elseif(to==2)
                    b2=bar(1:11,dumc,'FaceColor',lavender_scheme(colindx,:),'Facealpha',.1,'edgecolor','none','BarWidth', 1);
                    hold on
                else
                    b3=bar(1:11,dumi,'FaceColor',lime_scheme(colindx,:),'Facealpha',.1,'edgecolor','none','BarWidth', 1);
                    hold on
                end
    
            end
            ylim([0 1.1*max([dumn' dumc' dumi'])])
    
            pos_tmp = get(gca, 'Position');
            pos_list=[pos_list;pos_tmp];

            %if not top
            if(spix>no_treat_opts*no_ud_opts)
                   pos_list(end,2)=pos_list(spix-no_treat_opts*no_ud_opts,2)-pos_tmp(4);
            else
                %title(tolist(tl_idx),'FontSize',12);
                tl_idx=mod(tl_idx,3)+1;

            end

            %if 2nd
            if(mod(spix,3)==2)
                   pos_list(end,1)=pos_list(spix-1,1)+pos_tmp(3);
            end

            %if 3rd
            if(mod(spix,3)==0)
                %xpos ypos xwidth ywidth
                   pos_list(end,1)=pos_list(spix-2,1)+2*pos_tmp(3);
            end

            

            
            xlim([0.5 11.5])

            if(spix>no_treat_opts*no_ud_opts*(no_tps-1))
                xticks([1 6 11])
                xticklabels({'0' '0.5' '1'})
                %xlabel('x')
            else
                xticks([])
            end

            if(mod(spix,no_ud_opts*no_treat_opts)==1)
               % ylabel('Cell count')
            end

            if(mod(spix,no_ud_opts*no_treat_opts)==0)
               %yyaxis right
               %ylabel(daylist(dl_idx));
               %dl_idx=dl_idx+1;
               %yticks([])
            end

            if(mod(spix,3)==1)
                   
                   yticks([max([dumn' dumc' dumi'])/2 max([dumn' dumc' dumi'])]);
                   %yticks([max([dumn' dumc' dumi'])]);'
                   %yticks([max([dumn' dumc' dumi'])]);
                                     
            else
                   yticks([]);
            end

            


            set(gca, 'Position',[pos_list(end,1) pos_list(end,2) pos_tmp(3) pos_tmp(4)],'FontSize',12)

            set(findall(gcf,'type','text'), 'FontSize', fsz, 'Color', 'k','FontName', ft)
             set(gca,'FontSize', fsz, 'FontName', ft)

            spix=spix+1; 
        end

    end
end
