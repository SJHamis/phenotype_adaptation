% This code performs consistency analysis (CA)
% For a step-by-step walk-through of the method, see e.g.,
% https://doi.org/10.1142/9789811223495_0001

%Convention: statistical significance is
%small if A_hat in [0.5, 0.56)
%medium if A_hat in [0.56, 0.64)
%large ifA_hat in [0.64, 0.71)

load('consistency_analysis_workspace')
no_bins = 20;
no_reps = 100;
t_star = 28+28;
no_cells_matrix = zeros(no_treat_opts,no_ud_opts,no_bins*no_reps);
consistency_matrix = zeros(no_treat_opts,no_ud_opts,no_bins);
consistency_matrix_max = zeros(no_treat_opts,no_ud_opts);

% Compute the number of cells at time t_star for each simulation
for to = 1:no_treat_opts
    for uo = 1:no_ud_opts
        for r = 1:no_bins*no_reps
            for x_opt = 1:11
                tmp_x_conc = phenotype_distribution_matrix(to,uo,r,t_star,x_opt);
                no_cells_matrix(to,uo,r) = no_cells_matrix(to,uo,r) + tmp_x_conc;
            end
        end
    end
end

% Divide the simulations into no_bins = 20 bins
% Perform consistency analysis, read about it here:
for to = 1:no_treat_opts
    for uo = 1:no_ud_opts
        x0=reshape(no_cells_matrix(to,uo,1:100),1,[]);
        tmp=zeros(no_bins-1,1);
        for bin = 2:no_bins
            r0 = (bin-1)*100+1;
            rf = r0+99;
            xj=reshape(no_cells_matrix(to,uo,r0:rf),1,[]);
            [A, A_hat]=getA_measure(x0,xj);
            consistency_matrix(to,uo,bin)=A_hat;
            tmp(bin-1)=A_hat;
        end
        tmp
        consistency_matrix_max(to,uo)=max(tmp);
    end
end

consistency_matrix_max
%Statistical significance is small for 
no_ss_small = nnz(consistency_matrix_max < 0.56)
no_ss_med = nnz(consistency_matrix_max < 0.64) - no_ss_small
no_others = prod(size(consistency_matrix_max)) - no_ss_small - no_ss_med


