function initial_cell_list = setInitialCellList_custom(no_cells,no_data_cols,fitspace,dose_vec,uo,gr_x0,gr_x1,ic_factor)

    initial_cell_list = zeros(no_cells,no_data_cols);

    % Set cell index in column 1:
    initial_cell_list(:,1)=(1:no_cells)'; 

    % Second column. Parent is 0. Do nothing.

    % Populate phenotypic state in column 4:
    % For the directed strategies (3 and 4)
    if(uo>=3)
        initial_cell_list(:,4) = ic_factor*ones(no_cells,1);%Set x from a normal distibution.
    else
        mean_x_init= 0;
        std_x_init = 5;
        pdN = makedist('Normal',mean_x_init,std_x_init); %Halfnormal
        pdHN = truncate(pdN,0,10); 
        Y = random(pdHN,no_cells,1);      
        min_val = min(Y);
        max_val = max(Y);
        Y = (Y- min_val) / (max_val - min_val);
        Y = round(Y,1);
        initial_cell_list(:,4)=Y;
    end

    % Set phenotype state-index in column 3:
    initial_cell_list(:,3)=int64(10*initial_cell_list(:,4))+1;

    % Get fitness as fnc of drug and x
    for idx = 1:no_cells 
        initial_cell_list(idx,5)=getProliferativeFitness(0,initial_cell_list(idx,3),fitspace,dose_vec,gr_x0,gr_x1); 
    end

end