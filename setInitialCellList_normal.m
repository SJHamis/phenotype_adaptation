function initial_cell_list = setInitialCellList(no_cells,no_data_cols,fitspace,x0_opt,dose_vec)

    mean_x_init = 0.5;
    std_x_init = 0.1;

    initial_cell_list = zeros(no_cells,no_data_cols);

    % Set cell index in column 1:
    initial_cell_list(:,1)=(1:no_cells)'; 

    % Second column. Parent is 0. Do nothing.

    % Randomise phenotypic state in column 4:
    initial_cell_list(:,4) = normrnd(mean_x_init,std_x_init,no_cells,1)';%Set x from a normal distibution.
    initial_cell_list(:,4) = round(initial_cell_list(:,4)/(2*mean(initial_cell_list(:,4))),1);%Normalise x to be in [0,1].Then round to 1 dp.

    % Set phenotypic state-index in column 3:
    initial_cell_list(:,3)=int64(10*initial_cell_list(:,4))+1;

    % Get fitness as fnc of drug and x
    for idx = 1:no_cells 
        initial_cell_list(idx,5)=getProliferativeFitness(0,initial_cell_list(idx,3),fitspace,x0_opt,dose_vec); 
    end

end