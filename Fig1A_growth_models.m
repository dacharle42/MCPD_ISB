%Daniel Charlebois - Matlab R2017b - March 2018
%Matlab script for generating exponential, logistic, and Allee results
%for Figure 1. 

close all; clear all; clc;

%% Parameters
global r K Nc
N0 = 5*10^5; N02 = 1.4*10^7; %initial number of cells
r = 0.24;                    %population growth rate (per hour)
K = 10^7;                    %carrying capacity of the environment (logistic growth)
Nc = 4*10^5; Nc2 = 10^7;     %critical population size required for growth (Allee effect) 
dt = 0.1;                    %sampling interval (hours)
t_end = 48;                  %duration of experiment (hours)
t_array = 0:dt:t_end;        %time array (hours)

%% Exponential growth
%analytically solve equation
% syms N(t) g
% eqn = diff(N,t) == g*N;
% exp_eqn = dsolve(eqn);

%generate results for exponential analytic solution
i = 1; hld = size(t_array); N_exp = zeros(1,hld(2));
for t = t_array
   N_exp(i) = N0*exp(r*t);
   i = i + 1;
end

%% Logistic growth
% syms N(t) g K
% eqn = diff(N,t) == g*N-g*N^2/K;
% dsolve(eqn)

%generate results for logistic analytic solution
i = 1; hld = size(t_array); N_logis = zeros(1,hld(2)); N_logis2 = zeros(1,hld(2));
for t = t_array
   N_logis(i) = (K*N0)/(N0+(K-N0)*exp(-r*t));
   N_logis2(i) = (K*N02)/(N02+(K-N02)*exp(-r*t));
   i = i + 1;
end

%% Logistic growth
%Allee Effect
% syms N(t) Ncrit g K
% eqn = diff(N,t) == g*N*(1-N/K)*(N/Ncrit-1);
% dsolve(eqn)

%generate results for logistic analytic solution
%call ODE solver (ODE_solver.m)
N_Allee = N0;
[t, X] = Fig1A_ODE_solver(t_end,dt,N_Allee); 
N_Allee_data_above_Nc = X(:,1);

N_Allee = N02; Nc = Nc2;
[t, X] = Fig1A_ODE_solver(t_end,dt,N_Allee);
N_Allee_data_below_Nc_above_K = X(:,1);

N_Allee = N0; Nc = Nc2;
[t, X] = Fig1A_ODE_solver(t_end,dt,N_Allee);
N_Allee_data_below_Nc_below_K = X(:,1);

%% Figure
hold on
plot(t_array,N_exp,'k-')
plot(t_array,N_logis,'b-')
plot(t_array,N_logis2,'b--')
plot(t_array,N_Allee_data_above_Nc,'m-'); 
plot(t_array,N_Allee_data_below_Nc_above_K,'m--'); 
plot(t_array,N_Allee_data_below_Nc_below_K,'m-.');
K_data(1:hld(2))= K;
plot(t_array,K_data,'r--')
hold off
xlabel('time (hours)'); ylabel('number of cells')
legend('Exponential','Logistic (N0 < K)','Logistic (N0 > K)', ...
    'Allee (Nc < N0 < K)','Allee (K < N0 < Nc)','Allee (0 < N0 < Nc & K)','Carrying capacity (K)')
axis([0 t_end 0 N02])