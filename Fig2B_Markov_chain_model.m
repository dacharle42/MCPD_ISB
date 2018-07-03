%Daniel Charlebois - April 2018 - Matlab R2017b
%Discreate Markov chain model for mother-bud states obtained from 
%microfluidics-microscopy experiments shown in Figure 2.

close all; clc;

tic 

%% parameters
PM = 0.012; PB = 0.01; 
PS = 0.15; PD = 0.85; 
n = 10^3; 

%% initial values
S11 = 0.45; S10 = 0.16; S01 = 0.33; S00 = 0.06;
Si = [S11 S10 S01 S00]

%% Markov chain model
P = [(1-PM)*(1-PB) (1-PM)*PB PM*(1-PB) PM*PB; ...
     0 (1-PM) 0 PM; ...
     0 0 (1-PB) PB; ...
     PS 0 0 PD];
 
Sf = Si*P^n

%% figure
figure(1)
x = {'S_{11}' 'S_{10}' 'S_{01}' 'S_{00}'};
box on
set(gca,'fontsize',34)
h = bar(Sf); 
set(h,'linestyle','none')
set(gca,'XtickLabel',x(1,:))
ylabel('fraction of cells')
xlabel('phenotypic state')
xlabelHandle = get(gca, 'Xlabel');
set(xlabelHandle, 'FontSize', 11)

toc