%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% This MATLAB file produces the Proliferative Fitness Matrix (PFM) from 
% gowth rates that are extracted from in vitro data published in Kavran 
% et al (PNAS, 2022). 
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

close all
clear all

% Read in growth rate data from file.
growth_table = readtable('invitro_data/growth_rates.csv','VariableNamingRule','preserve');


% Make one table per drug-exposure setting:
% No drug exposure 
% 1 week  (500 nM) drug exposure
% 2 weeks (500 nM) drug exposure
% 3 weeks (500 nM) drug exposure
% old: Columns: 1 (previous drug exposure). 2 (rep). 3 (dose). 6 (growth rate).
% new: 1,2,5 9
column_filter = [1,2,5,9];
tab_v600E=growth_table(strcmp(growth_table.sample,'V600E'),column_filter);
tab_1wks=growth_table(strcmp(growth_table.sample,'Week 1'),column_filter);
tab_2wks=growth_table(strcmp(growth_table.sample,'Week 2 C'),column_filter);
tab_3wks=growth_table(strcmp(growth_table.sample,'Week 3 C'),column_filter);
tab_4wks=growth_table(strcmp(growth_table.sample,'Week 4 C'),column_filter);

no_reps = 4;
no_doses = size(tab_v600E,1)/no_reps;
dose_vec = unique(table2array(tab_v600E(:,3)));
%dose_vec(1)=0;

% Create and populate vectors with average growth rate data (over four 
% replicates). 1 vector per drug-exposure setting.
pf_v600E = zeros(no_doses,1);
pf_1wks = zeros(no_doses,1);
pf_2wks = zeros(no_doses,1);
pf_3wks = zeros(no_doses,1);
pf_4wks = zeros(no_doses,1);

for d = 1:no_doses
    start_row=(d-1)*no_reps+1;
    end_row=d*no_reps;
    pf_v600E(d) = mean(table2array(tab_v600E(start_row:end_row,4)));
    pf_1wks(d)  = mean(table2array(tab_1wks(start_row:end_row,4)));
    pf_2wks(d)  = mean(table2array(tab_2wks(start_row:end_row,4)));
    pf_3wks(d)  = mean(table2array(tab_3wks(start_row:end_row,4)));
    pf_4wks(d)  = mean(table2array(tab_4wks(start_row:end_row,4)));
end

% Combine mean vectors to a matrix.
pf_vecs=[pf_v600E pf_1wks pf_2wks pf_3wks pf_4wks];

figure
tl=tiledlayout(2,2);

%%%%%%%%%%%%%% SUBPLOT 1 %%%%%%%%%%%%%%
ft = 'Times';
fsz = 20;        

set(gca, 'ColorOrderIndex', 1);
% Plot mean growth rate values.
nexttile
set(gca, 'ColorOrderIndex', 1);
for d=1:size(pf_vecs,2)
    semilogx(dose_vec,pf_vecs(:,d),'LineWidth',2)
    hold on
end

% Pick out all (dose,growth rate) data pairs in order to scatter plot them.
s0=table2array(tab_v600E(:,3:4));
sw1=table2array(tab_1wks(:,3:4));
sw2=table2array(tab_2wks(:,3:4));
sw3=table2array(tab_3wks(:,3:4));
sw4=table2array(tab_4wks(:,3:4));

% Scatter-plot all growth rate values.
set(gca, 'ColorOrderIndex', 1);
hold on
scatter((s0(:,1)),s0(:,2),"^",'LineWidth',2)
hold on
scatter((sw1(:,1)),sw1(:,2),'LineWidth',2)
hold on
scatter((sw2(:,1)),sw2(:,2),'LineWidth',2)
hold on
scatter((sw3(:,1)),sw3(:,2),'LineWidth',2)
hold on
scatter((sw4(:,1)),sw4(:,2),'LineWidth',2)
hold on
semilogx(dose_vec,0*dose_vec,'k--','LineWidth',1.5)

%0 wks is V600E in Kavran et al. (2022)
legend('0 wks mean','1 wks mean','2 wks mean',...
    '3 wks mean','4 wks mean',...
    '0 wks','1 wks','2 wks','3 wks','4 wks','Location','southwest','FontSize',20);
xlabel('BRAFV600E-inhibitor dose (nM)')
ylabel('Growth rate (day^{-1})')
set(findall(gcf,'type','text'), 'FontSize', fsz, 'Color', 'k','FontName', ft)
set(gca,'FontSize', fsz, 'FontName', ft)
ylim([-0.9 0.3])
yticks(-0.8:0.2:0.2);
labelvec = dose_vec;
labelvec(1) = 0;
xticks(dose_vec)
xticklabels(labelvec)

%%%%%%%%%%%%%% SUBPLOT 2 %%%%%%%%%%%%%%
% Plot growth rate values for sensitive-to-resistant interpolated 
% phenotypes.
nexttile

% Define customised colour scheme to plot from blue to red.
mycol_blue2red= ...
  [0 0 255
  26 0 230
  51 0 204
  77 0 179
  102 0 153
  128 0 128
  153 0 102
  179 0 77
  204 0 51
  230 0 26
  255 0 0]/255;
colororder(mycol_blue2red);

% Get mean value for resistant (drug-dependent) phenotypes (those which have been exposed to
% drugs for 1-4 weeks.
mean_resistant=mean(pf_vecs(:,2:5),2);
semilogx(dose_vec,pf_vecs(:,1),'blue','LineWidth',2)
hold on
semilogx(dose_vec,mean_resistant,'red','LineWidth',2)
hold on
semilogx(dose_vec,0*dose_vec,'k--','LineWidth',1.5)
hold on


% Set the chosen number of phenotypes.
no_x = 11;

% Create and allocate values to the proliferative fitness matrix via linear
% interpolation. 
%format long
proffit_matrix = zeros(no_doses,no_x);
for d = 1:no_doses
    s=pf_vecs(d,1);
    r=mean_resistant(d);
    proffit_matrix(d,:)=linspace(s,r,no_x);
end

set(gca, 'ColorOrderIndex', 1);
colororder(mycol_blue2red);
for d=1:no_doses
    for x = 1:no_x
        hold on
        sized=10*d;
        scatter(dose_vec(d),proffit_matrix(d,x),sized,'filled')
    end
end


xlabel('BRAFV600E-inhibitor dose (nM)')
ylabel('Growth rate (day^{-1})')



% 500 nM interpolation
% Interpolate drug naive growth rate at 500 nM
x0 = 300;
delta_x = log(1000/300);
y0 = proffit_matrix(9,1);
delta_y = proffit_matrix(10,1)-y0;
slope = delta_y/delta_x;
gr_500_x0 = y0 + slope * log(500/x0);
% Interpolate drug resistant growth rate at 500 nM
y0 = proffit_matrix(9,end);
delta_y = proffit_matrix(10,end)-y0;
slope = delta_y/delta_x;
gr_500_x1 = y0 + slope * log(500/x0);


set(gca, 'ColorOrderIndex', 1);
colororder(mycol_blue2red);
vec_500nm = linspace(gr_500_x0,gr_500_x1,no_x);
for x = 1:no_x
    hold on
    sized=10*(9+10)/2;
    scatter(500,vec_500nm(x),sized,'o','LineWidth',2)
end

legend('Drug-naive phenotype','Drug-addicted phenotype','','','Location','southwest','FontSize',20);
set(findall(gcf,'type','text'), 'FontSize', fsz, 'Color', 'k','FontName', ft)
set(gca,'FontSize', fsz, 'FontName', ft)
ylim([-0.9 0.3])
yticks(-0.8:0.2:0.2);
dose_vec_mod = [dose_vec(1:9);500;dose_vec(10:end)];
labelvec_mod = [labelvec(1:9);500;labelvec(10:end)];
xticks(dose_vec_mod)
xticklabels(labelvec_mod)

%%%%%%%%%%%%%% SUBPLOT 3 %%%%%%%%%%%%%%
% Plot the phenotype and drug values that map to the proliferative fitness 
% matrix.
nexttile
semilogy(0,0)
xlim([0 1])
ylim([0 10000])
hold on
for d=1:no_doses
    set(gca, 'ColorOrderIndex', 1);
    colororder(mycol_blue2red);
    for x=1:no_x
        sized=10*d;
        x_val=(x-1)*0.1;
        scatter(x_val,dose_vec(d),sized,'filled')
    end
end
ylabel('BRAFV600E-inhibitor dose (nM)')
xlabel('Phenotypic state value (x)')

set(findall(gcf,'type','text'), 'FontSize', fsz, 'Color', 'k','FontName', ft)
set(gca,'FontSize', fsz, 'FontName', ft)

yticks(dose_vec)
yticklabels(labelvec)

%%%%%%%%%%%%%% SUBPLOT 4 %%%%%%%%%%%%%%
% Plot values of the proliferative fitness matrix in a heatmap.
nexttile
proffit_matrix_ud=flipud(proffit_matrix);
xvalues = string(linspace(0,1,no_x));

yvalues = flipud(num2str(labelvec));

hold off
h=heatmap(xvalues,yvalues,proffit_matrix_ud,'Colormap',gray);
ylabel('BRAFV600E-inhibitor dose (nM)')
xlabel('Phenotypic state value (x)')



set(findall(gcf,'type','text'), 'FontSize', fsz, 'Color', 'k','FontName', ft)
set(gca,'FontSize', fsz, 'FontName', ft)


%fontsize(tl,16,"points")

% Save the PFM and vector of drug doses
writematrix(proffit_matrix,'PFM.xls')
writematrix(dose_vec,'dosevector.xls')