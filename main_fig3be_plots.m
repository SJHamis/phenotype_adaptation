
fsz=10;

for opt = 1:2
    if(opt==1)
       avg_rates=load('avg_rates_56h_ic0.mat');
    else
       avg_rates=load('avg_rates_56h_ic1.mat');
    end
    avg_rates=cell2mat(struct2cell(avg_rates));

    daily_pheno_ups_list = 1:10;%1:10;
    delta_treat_list =1:10;% 1:11%[1 2 7 14 28];
    max_delta_treat = length(delta_treat_list);
    
    g1 = fitness_matrix(1,1);
    g3 = fitness_matrix(1,11);

for sp = 1:2

    % Prepare figure for svg saving.
    fig = figure; 
    fig.Renderer = 'painter';
    fig.Position = [0 0 300 300];
    tl=tiledlayout(2,2);
    tl.TileSpacing = 'loose';
    ft = 'Arial';
    fsz = 18; 

    tmp_map = avg_rates(sp,:,:);
    tmp_map=reshape(tmp_map,length(daily_pheno_ups_list),max_delta_treat);
    data=tmp_map;

    if(sp==1)%300 nM
        g2 = fitness_matrix(9,1);
        g4 = fitness_matrix(9,11);
    else %500 nM
        g2=-0.2866; %interpolated
        g4=0.0741;
    end

   lambda_target= g4;
    delta_h=0.1;   
    Tw_more = (-g1+g2+g3-g4) / (4 * lambda_target - 2*(g1+g4));
    nw_more = Tw_more/delta_h;
    eff_gr_cont = g4;

    if(sp==1)
        mind=-0.03;
    else
        mind=-0.08;
    end
    
    n_colors = 256; % Number of colors in the colormap
    % Specify the values corresponding to these colors
    if(eff_gr_cont< max(data(:)))
        cmap = [255,223,0 %yellow
        247,247,247 %white
        0,0,0]/255; %black
        max_gr = 0.096; 
        max(data(:))
        color_values = [mind, eff_gr_cont, max_gr];
    else
        min(data(:))
        cmap = [255,223,0 %yellow
            247,247,247]/255; %white
        color_values = [mind, eff_gr_cont];
    end

    % Interpolate the colormap for smooth transitions

    if(eff_gr_cont< max(data(:)))
    new_cmap = interp1(color_values, cmap, linspace(min(color_values), max(color_values), n_colors));
    else
    new_cmap = interp1(color_values, cmap, linspace(min(color_values), eff_gr_cont, n_colors));
    end

    
    % Create the heatmap
    imagesc(data);
    axis square;
    colorbar; % Display colorbar
    
    % Apply the custom colormap
    colormap(new_cmap);
    clim([min(color_values), max(color_values)]); % Ensure the color limits match the colormap
    
    % Add a colorbar
    cb = colorbar;
    set(gca, 'FontSize', 14);


    if(sp==1)
        lab_values = [-0.02 0 0.02 0.05 eff_gr_cont];
    else
        lab_values = [-0.05 0 0.01 0.05 eff_gr_cont];
    end
    
    % Specify positions for custom labels (match these with data or color range)
    tick_positions = lab_values;  % Positions on the colorbar corresponding to color_values

    %tick_labels = {'lambda_gr'};  % Custom labels for these positions
    
    % Set the colorbar ticks and labels
    set(cb, 'Ticks', tick_positions,'FontSize',12)%, 'TickLabels', tick_labels);
    
    % Optional: Labeling axes and adding title
    xlabel('interval length in days ({\itT})','FontSize',14);
    ylabel('phenotype updates per day ({\itη})','FontSize',14);
    
    % Hold the plot to overlay the curve
    hold on;
    
    % Define the range for X-axis and calculate Y based on XY = 33
    X = linspace(min(get(gca, 'XLim')), max(get(gca, 'XLim')), 1000); % X values based on current axes limits
    x= 0:100;
    y=nw_more ./ x; 
    Y = nw_more ./ X; % Calculate Y for XY = 33
    
    % Plot the curve on top of the heatmap
    plot(X, Y, 'g--', 'LineWidth', 3); % 'r' specifies the color red and 'LineWidth' sets the line width
    set(gca,'ColorScale','log')

    title=strcat('Figure_3_xICopt',num2str(opt),'_dose',num2str(sp),'.svg');
    saveas(fig, title); 
end

end
