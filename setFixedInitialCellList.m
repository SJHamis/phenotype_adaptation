function initial_cell_list = setFixedInitialCellList(no_cells,no_data_cols,fitspace,x0_opt,dose_vec)

    initial_cell_list = zeros(no_cells,no_data_cols);
    % Set cell index in column 1:
    initial_cell_list(:,1)=(1:no_cells)'; 
   
    % Set ancestor to be itself in column 2:
    initial_cell_list(:,2)=(1:no_cells)';

    % Fix phenotypic state in column 4:
    initial_cell_list(:,4) = linspace(0,1,no_cells)';%Fix x

    % Fix phenotypic state-index in column 3:
    initial_cell_list(:,3)=int64(10*initial_cell_list(:,4))+1;

    % Get fitness as fnc of drug and x
    for idx = 1:no_cells 
        initial_cell_list(idx,5)=getProliferativeFitness(0,initial_cell_list(idx,3),fitspace,x0_opt,dose_vec); 
    end

end