addpath perfusion_classes;
clear all;
clc;

N=100000;
seed=1;
solver = PerfusionSolver(N,seed);

y_pos = linspace(20,370,N); %Linear spacing between 20 - 370 micrometers
x_pos = rand(1,N).*440.*ones(1,N)+20; %Randomly distribute them from 20 - 460 micrometers
z_pos = 0*ones(1,N); %All particles at z=20 micrometers

solver.setInitPos([x_pos',y_pos',z_pos']);

%Set Gradient Profile
dt = linspace(0,1,1000);
gG = -1 * ones(1,1000);
gG(dt>=500)=1;
T = 1000;
Gmax = 1;

solver.setGradient(struct('dt',dt,'gG',gG),Gmax,T);
solver.setLengthDistribution("Constant",1);
solver.setAnglesDistributions("Isotropic");

[phase,capillaries] = solver.computePerfusion("parallel");