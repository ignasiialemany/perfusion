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
[normalized_phase,~] = solver.computePerfusion("serial");
[phase,b_value,G_max] = solver.computePhase(normalized_phase,"Gmax",35.690e-6);
tol=0.001;
%Phase is known for bipolar gradient = -gamma*G_max*v*T^2*0.25
assert(max(phase(:,3)+954.786)<tol);


%% Test 2 - Check bipolar gradient for isotropic directions comparing it to 
% Flow-Compensated Intravoxel Incoherent Motion
% Diffusion Imaging - Andreas Wetscherek,1* Bram Stieltjes,2 and Frederik Bernd Laun1,2

N_segment=1;
solver.setLengthDistribution("Constant",N_segment);
solver.setAnglesDistributions("Isotropic");
[normalized_phase,~] = solver.computePerfusion("serial");
comp_values = load("normalized_phase_test2_N1.mat");
assert(isequal(normalized_phase,comp_values.normalized_phase));

N_segment=3;
solver.setLengthDistribution("Constant",N_segment);
solver.setAnglesDistributions("Isotropic");
[normalized_phase,capillaries] = solver.computePerfusion("serial");
comp_values = load("normalized_phase_test2_N3.mat");
assert(isequal(normalized_phase,comp_values.normalized_phase));






