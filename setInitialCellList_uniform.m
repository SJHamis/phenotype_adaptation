function initial_cell_list = setInitialCellList_uniform(no_cells,no_data_cols,fitspace,dose_vec)

    initial_cell_list = zeros(no_cells,no_data_cols);

    % Set cell index in column 1:
    initial_cell_list(:,1)=(1:no_cells)'; 

    % Second column. Parent is 0. Do nothing.

    % Fix phenotypic state in column 4:
    for idx = 1:no_cells 
        initial_cell_list(idx,4) = mod((idx-1)*0.1,1.1);
    end
    
    % Set phenotypic state-index in column 3:
    initial_cell_list(:,3)=int64(10*initial_cell_list(:,4))+1;

    % Get fitness as fnc of drug and x
    for idx = 1:no_cells 
        initial_cell_list(idx,5)=getProliferativeFitness(0,initial_cell_list(idx,3),fitspace,dose_vec); 
    end

end