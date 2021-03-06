% Test Perfusion Code
addpath ../perfusion_classes;
%Common Parameters
N=10000;
seed=1;

%% Common Part
solver = PerfusionSolver(N,seed);
dt = linspace(0,1,1000);
gG = -1 * ones(1,1000);
gG(dt>=0.5)=1;
T = 1000;
solver.setGradientProfile(struct('dt',dt,'gG',gG),T);

%% Test 1 - Check bipolar gradient for 0 zenith angle (only z-axis movement)
solver.setCapillaryDistributions("Constant","Isotropic",1);
solver.setVelocityValue(0.4);
% %Constant Zenith Angle
solver.direction_func = @(stream,mean) [0,0,1];
[normalized_phase,~] = solver.computePerfusion("serial");
[phase,b_value,G_max] = solver.computePhase(normalized_phase,"Gmax",35.690e-6);
tol=0.001;
% %Phase is known for bipolar gradient = -gamma*G_max*v*T^2*0.25
assert(max(phase(:,3)+954.786)<tol);


%% Test 2 - Check bipolar gradient for isotropic directions through paper comparison
% Flow-Compensated Intravoxel Incoherent Motion
% Diffusion Imaging - Andreas Wetscherek,1* Bram Stieltjes,2 and Frederik Bernd Laun1,2
N_segment=1;
solver.setCapillaryDistributions("Constant","Isotropic",N_segment);
[normalized_phase,~] = solver.computePerfusion("parallel");

N_segment=3;
solver.setCapillaryDistributions("Constant","Isotropic",N_segment);
[normalized_phase,capillaries] = solver.computePerfusion("parallel");




