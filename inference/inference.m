%% LOAD SYNTHETIC DATA (WILL EVENTUALLY LOOP OVER PUP'S)

Ycont = readmatrix('../output_data_cont.xls');
Yint = readmatrix('../output_data_cont.xls');

% Join [Ycont,Yint], each row is a rep (just three rows for now, say)
Ydata = [Ycont,Yint]; Ydata = Ydata(1:3,:);

%% SETUP LIKELIHOOD FUNCTION AND LOG POSTERIOR
% - Assume normal on the log scale
% - Uniform priors, up to a constant

% Specify model, load prior
model = @ode1; ode1_prior;
npars = length(lower);

% Log posterior
llpost = @(expars) logpost(expars,model,Ydata,lower,upper);

%% INFERENCE USING ADAPTIVE MCMC

% Initial guess
p0 = [1.0,lower + rand() * (upper - lower)];

% Todo: find a good adapative mcmc package in matlab