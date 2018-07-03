%Daniel Charlebois - Winter 2018 - Matlab R2017b
%This script generates growth curves for a deterministic and a
%stochastic version of the Baranyi model, as shown in Figure 1B.
%There is also an ODE solver version (commented by default) of the Baranyi 
%model included.

clc; close all;

tic

%% parameters
global mu N_max lambda n
mu = 0.24;         %growth rate (per hour)
N_max = 10^6;      %carrying capacity
lambda = 4;        %lag time (hour)
n = 4;             %parameters for adaptation function
n_runs = 10;       %number of runs
sigma_0 = 0;       %no noise 
sigma = 0.035;     %SDE noise strength
n0 = 10^4;         %initial number of cells
t_end = 48;        %simulation time (hours)
n_step = 96;       %number of steps                
h = t_end/n_step;  %time step
t=(0:h:t_end);     %time array

%% SDE solver
size_t = size(t);
n_data_0 = zeros(1,size_t(2));
n_data = zeros(n_runs,size_t(2));
n_data_0(1,1) = n0;
for j = 1:n_runs
    
    n_data(j,1) = n0;
    
    for i=1:n_step
        if j == 1
            n_data_0(j,i+1) = n_data_0(j,i) + (mu*n_data_0(j,i)*(t(i)^n/(lambda^n+t(i)^n)))*(1-(n_data_0(j,i)/N_max))*h + sigma_0*n_data_0(j,i)*randn*sqrt(h);
        end
        n_data(j,i+1) = n_data(j,i) + (mu*n_data(j,i)*(t(i)^n/(lambda^n+t(i)^n)))*(1-(n_data(j,i)/N_max))*h + sigma*n_data(j,i)*randn*sqrt(h);
    end
    
end

% dt=1; n_ode=n0;
% %% ODE solver
% [t_ODE, X] = Fig1B_ODE_Baranyi(t_end,dt,n_ode);

%% save data
save('baranyni.mat','t','n_data_0','n_data');
 
%% plot
figure; 
hold on
plot(t,n_data_0,'r-','LineWidth',4);
plot(t,n_data,'b--','LineWidth',2);
% plot(t_ODE,X,'g-','LineWidth',2)
hold off
xlabel('time (hours)'); ylabel('number of cells');
legend('Baranyi noise strength = 0','Baranyi noise strength = 0.035')

toc