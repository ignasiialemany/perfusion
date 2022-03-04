clear all
clc

addpath external/yamlmatlab
addpath perfusion_classes/

%Perfusion simulation

[dt,gG,g_strength,T]  = create_seq("sequence.yml");

%gG = (mT/micrometer)
%T = (m/s)

N_spins=100000;

%We seed the spins in a smaller voxel so that we make sure almost all of
%them stay within the scanning voxel (3x3x8 mm^3)

space_traveled = 0.5*T; %(v=0.5 micrometers/ms) + safety factor of 1.1

k=[3.25];
for k_val = k
    seed=floor(unifrnd(1,10000));
    
    %Check that seeding z=0 and seeding z to all the voxel is the same
    %X-Y plane 
    x = unifrnd(-space_traveled,3000+space_traveled,N_spins,1);
    y = unifrnd(-space_traveled,3000+space_traveled,N_spins,1);
    z = unifrnd(-space_traveled,8000+space_traveled,N_spins,1);
    
    %x = unifrnd(-space_traveled,3000+space_traveled,N_spins,1);
    %y = unifrnd(-space_traveled,3000+space_traveled,N_spins,1);
    %z = unifrnd(-space_traveled,8000+space_traveled,N_spins,1);
    %z = unifrnd(0,8000,N_spins,1);
    solver = PerfusionSolver(N_spins,seed);
    %Init position and capillary distributions
    init_pos = [x,y,z];
    solver.setInitPos(init_pos);
    solver.setGradientProfile(struct('dt',dt,'gG',gG),T);  
    solver.setCapillaryDistributions("Weibull","Anisotropic-VonMises",k_val);
    solver.isFibreDirection=true;
    name = "FibreDirection_VonMises_" + num2str(k_val);
    [nP,cp,final_pos] = solver.computePerfusion("parallel");
    save(name);   
end








