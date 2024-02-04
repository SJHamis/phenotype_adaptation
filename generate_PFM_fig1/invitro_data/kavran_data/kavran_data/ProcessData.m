%% PROCESSDATA
%  Download and then process data to obtain the data in equivalent growth
%  rate form. 

%% DOWNLOAD DATA
DownloadData;

%% OBTAIN CONTINUOUS THERAPY GROWTH RATE FOR NORMALISATION

% Read data
data_gr = readtable("kavran_growthcurve.csv",ReadVariableNames=true);

% Subset to obtain V600E under continuous treatment
idx = strcmp(data_gr.sample,'V600E') .* strcmp(data_gr.schedule,'Continuous');

% Convert to matrix, and then a vector of cell counts (day 7 onwards)
cell_counts_mat = table2array(data_gr(idx == 1,5:end));
cell_counts_vec = cell_counts_mat(:);

% Day number to a vector
day_mat = repmat(7:7:28,6,1);
day_vec = day_mat(:);

% Exclude the outlier
idx = cell_counts_vec < 7000;
cell_counts_vec = cell_counts_vec(idx);
cell_counts_mat(cell_counts_mat > 7000) = NaN;
day_vec = day_vec(idx);

% Perform regression to get the growth rate

    % Design matrix
    X = [ones(size(day_vec)) day_vec];

    % Response
    Y = log(cell_counts_vec);

    % Estimate for the slope
    beta = X \ Y;
    lambda_cont = beta(2);

    % Approximate confidence intervals
    s2 = sum((X * beta - Y).^2) / (length(Y) - 2);
    beta_sig = s2 * inv(X' * X);
    lambda_ci = lambda_cont + norminv(0.975) * sqrt(beta_sig(2,2)) * [-1,1];
    
% Produce supplementary plot
clf;

% Data (error bars)
errorbar(7:7:28,nanmean(cell_counts_mat),nanstd(cell_counts_mat), ...
    '.','MarkerSize',40,'LineWidth',2,'Color',[0.8,0.5,0.0]); hold on;

% Data (raw data)
scatter(day_vec,cell_counts_vec,50,'d','filled','LineWidth',1.5,'Color',[0.7,0.2,0.0]); hold on;

% Fit
fplot(@(t) exp(beta(1) + t * beta(2)),'k--','LineWidth',2)
xlim([0.0,30.0]); xticks(0:7:28); grid on;
xlabel('Time [d]'); ylabel('Cell Count');

% Save
saveas(gcf,'Fig_S1_V600E_Continuous_Fit.png');

%% DOSE RESPONSE TO GROWTH RATES

% Read data, convert to long form
data_dr_raw = readtable("kavran_doseresponse.csv",ReadVariableNames=true);
[m,n] = size(data_dr_raw); n = n - 2;

    % Obtain doses
    dose_str = string(data_dr_raw.Properties.VariableNames(3:end));
    dose_dbl = str2double(strrep(strrep(dose_str,"_","."),"x",""));
    dose_mat = repmat(dose_dbl,m,1);
    dose = dose_mat(:);
    log_dose = log(dose);

    % Obtain measurements (fc)
    fc = table2array(data_dr_raw(:,3:end));
    fc = fc(:);

    % Obtain Sample names
    sample = repmat(string(data_dr_raw.sample),1,n);
    sample = sample(:);

    % Obtain reps
    rep = str2double(repmat(string(data_dr_raw.rep_number),1,n));
    rep = rep(:);

    % New variables: treatment, and day
    sample_vals = unique(sample);
    day_repl_vals = [0,0,7,14,14,21,21,28,28];
    trt_rep_vals = {'nil','nil','continuous','continuous','intermittent',...
        'continuous','intermittent','continuous','intermittent'};
    day = zeros(m*n,1);
    treatment = sample;
    for i = 1:length(sample_vals)
        idx = sample == sample_vals(i);
        day(idx) = day_repl_vals(i);
        treatment(idx) = trt_rep_vals{i};
    end

% Unnormalised growth rates
growth_rate_unnormalised = log(fc) / 3;

% Obtain (unnormalised) mean for continuously treated at 500nM dose

    % Mean for 300nM
    unorm_lambda_cont_300nM = mean(growth_rate_unnormalised(...
        ((dose == 300) .* strcmp(treatment,'continuous')) == 1));

    % Mean for 1000nM
    unorm_lambda_cont_1000nM = mean(growth_rate_unnormalised(...
        ((dose == 1000) .* strcmp(treatment,'continuous')) == 1));

    % Mean for 500nM
    unorm_lambda_cont_500nM = unorm_lambda_cont_300nM + (unorm_lambda_cont_1000nM - unorm_lambda_cont_300nM) * (500 - 300) / (1000 - 300);

% Normalised growth rates
growth_rate = growth_rate_unnormalised - unorm_lambda_cont_500nM + lambda_cont;

% Construct new table
data_dr = table(sample,rep,day,treatment,dose,log_dose,fc,growth_rate_unnormalised,growth_rate);

% Remove data for vector and save table
data_dr = data_dr(~strcmp(data_dr.sample,"Vector Control"),:);
writetable(data_dr,"growth_rates.csv");

% Produce plot
clf;
log_dose_vals = unique(data_dr.log_dose);
gr_naive = zeros(size(log_dose_vals));
gr_addict = zeros(size(log_dose_vals));
fc_naive = zeros(size(log_dose_vals));
fc_addict = zeros(size(log_dose_vals));
for i = 1:length(log_dose_vals)
    gr_naive(i) = mean(data_dr.growth_rate(...
        ((data_dr.log_dose == log_dose_vals(i)) .* strcmp(data_dr.treatment,'nil')) == 1));
    gr_addict(i) = mean(data_dr.growth_rate(...
        ((data_dr.log_dose == log_dose_vals(i)) .* strcmp(data_dr.treatment,'continuous')) == 1));
    fc_naive(i) = mean(data_dr.fc(...
        ((data_dr.log_dose == log_dose_vals(i)) .* strcmp(data_dr.treatment,'nil')) == 1));
    fc_addict(i) = mean(data_dr.fc(...
        ((data_dr.log_dose == log_dose_vals(i)) .* strcmp(data_dr.treatment,'continuous')) == 1));
end

subplot(1,2,1);
semilogx(exp(log_dose_vals),fc_naive,'b','LineWidth',2); hold on;
semilogx(exp(log_dose_vals),fc_addict,'r','LineWidth',2);
plot([500,500],[0.0,1.0],'k:','LineWidth',2);
legend('Naive','Addict','Location','SouthWest');
xlabel('log dose [nM]'); ylabel('fc')

subplot(1,2,2)
semilogx(exp(log_dose_vals),gr_naive,'b','LineWidth',2); hold on;
semilogx(exp(log_dose_vals),gr_addict,'r','LineWidth',2)
semilogx([min(dose),max(dose)],[0.0,0.0],'k--','LineWidth',2)
plot([500,500],[-1.0,0.3],'k:','LineWidth',2);
ylim([-1.0,0.3]);
legend('Naive','Addict','Location','SouthWest')
xlabel('log dose [nM]'); ylabel('fc')

% Save
saveas(gcf,'Fig_S2_FC_to_GrowthRates.png');