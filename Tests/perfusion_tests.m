% Test Perfusion Code
addpath perfusion_classes;
clear all;
clc;

%Common Parameters
N=100;
seed=1;

%% Test 1 - Check bipolar gradient for 0 zenith angle (only z-axis movement)
solver = PerfusionSolver(N,seed);
dt = linspace(0,1,1000);
gG = -1 * ones(1,1000);
gG(dt>=0.5)=1;
T = 1000;
solver.setGradientProfile(struct('dt',dt,'gG',gG),T);
solver.setLengthDistribution("Constant",1);

%Constant Zenith Angle
solver.zenith_func = @(x,mean) 0;
solver.azimuth_func = @(x) 2*pi*x;
[normalized_phase,capillaries] = solver.computePerfusion("serial");
[phase,b_value,G_max] = solver.computePhase(normalized_phase,"Gmax",35.690e-6);
tol=0.001;
%Phase is known for bipolar gradient = -gamma*G_max*v*T^2*0.25
assert(max(phase(:,3)+954.786)<tol);





