%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% This MATLAB file checks which doses are more 
% effective when given in intermittent schedules 
% and estimates proliferation rates for 500 nM doses.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% clear all
% close all

PFM1=readmatrix('generate_PFM_fig1/PFM.xls');
dose_vec=readmatrix('generate_PFM_fig1/dosevector.xls');

%In the PFM,
% row 9: dose 300.
% row 10: dose 1000.
% Use linear interpolation to estimate growth rates for 500 nM

x0 = 300;
delta_x = log(1000/300);

% Interpolate drug naive growth rate at 500 nM
y0 = PFM1(9,1);
delta_y = PFM1(10,1)-y0;
slope = delta_y/delta_x;
gr_500_x0_raw = y0 + slope * log(500/x0)

% Interpolate drug resistant growth rate at 500 nM
y0 = PFM1(9,end);
delta_y = PFM1(10,end)-y0;
slope = delta_y/delta_x;
gr_500_x1_raw = y0 + slope * log(500/x0)
