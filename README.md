# phenotype_adaptation
Code base for the Phenotype Adaptation project.

Run MATLAB files starting with the key-word “main_” to produce manuscript result figures. 

1. The MATLAB file generate_PFM_fig1/main_generate_PFM_fig1.m generates the Proliferative Fitness Matrix (which gets saved as PFM.xls) and Figure 1.

2.  The plots in Figures 2b and 2c are, respectively, generated with the MATLAB files main_fig2b_make_example_phenoupdates.m and main_fig2c_phenotype_bar_dynamics.m.

3. Figure 3 plots are main_fig3_ContVSInt.m. Mean cell count data from the simulation experiments are saved in “treatment_responses_cont.xls” (for continous treatment) and “treatment_responses_int.xls” for intermittent treatment.

In these xls files, each column corresponds to a time point (day 1,2,..,28).
Each row (1,2,...20) corresponds to a drug dose (30, 100, 300, 500, 1000) nM and phenotype update (PUP1,2,3,4) model. 
So row R in the xls files has 
dose-index ceil(R/4), 
PUP-index mod(R-1,4)+1.

xls file structure:
Row, dose-index, PUP-index
1 1 1
2 1 2
3 1 3
4 1 4
5 2 1 
6 2 2
7 2 3
8 2 4
9 3 1
...
19 5 3
20 5 4
