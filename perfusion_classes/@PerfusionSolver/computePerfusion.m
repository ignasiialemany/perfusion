function [phase,capillaries] = computePerfusion(obj,execution_mode,varargin)

%Inputs
inputs = obj.inputs;
num_streams = inputs.n_particles;
seed = inputs.seed;
pos = inputs.init_pos;
gG = inputs.gradient.gG;
dt = inputs.gradient.dt;
T = inputs.T;
v = inputs.velocity_value;

%Output
phase = zeros(size(pos));
capillary_info=struct();
capillary_info.direction = [];
capillaries(1,inputs.n_particles)=capillary_info();

%Parallel code, create workers if not initialized.
if strcmp(execution_mode,"parallel")
    poolobj = gcp('nocreate');
    %Retrieve Number of Cores
    if numel(varargin)>0
        n_cores = varargin{1};
    else
        n_cores = maxNumCompThreads;
    end
    %Setup pool of workers if it is already not created
    if isempty(poolobj)
        parpool(n_cores);
    end
    
    parfor i=1:inputs.n_particles
        %Create stream for each particle
        stream = RandStream.create('mlfg6331_64', 'Seed', seed, ...
            'NumStreams', num_streams, 'StreamIndices', i);
        pos_particle = pos(i,:);
        [phase(i,:),capillary_directions] = phase_one_particle(obj,stream,pos_particle,T,v,dt,gG);
        capillaries(1,i).direction = capillary_directions;
    end
else
    %Loop
    for i=1:inputs.n_particles
        %Create stream for each particle
        stream = RandStream.create('mlfg6331_64', 'Seed', seed, ...
            'NumStreams', num_streams, 'StreamIndices', i);
        pos_particle = pos(i,:);
        [phase(i,:),capillary_directions] = phase_one_particle(obj,stream,pos_particle,T,v,dt,gG);
        capillaries(1,i).direction = capillary_directions;
    end
    
end
phase = phase .* (1/inputs.gradient.a);
end

function [phase,capillary_directions] = phase_one_particle(obj,stream,pos,T,v,dt,gG)
t_init=0;
t_end=0;
capillary_directions=[];
phase=zeros(1,3);
while(t_end<T)
    
    %Compute capillary segment
    length_value = obj.length_func(stream);
    
    %Compute zenith angle according to y_spacing / Virtual tissue
    %recreation
    y_spacing = 392.04; %y-spacing of block 
    block_index = floor(pos(2)/y_spacing);
    angle = pi/2 - deg2rad(3.9204)*block_index;
    
    %Get the mean fibre direction. Azimuth angle is not important is set to
    %0
    default_azimuth=0.;
    [dir1,dir2,dir3] = sph2cart(deg2rad(default_azimuth),angle,1);
    direction = obj.direction_func(stream,[dir1,dir2,dir3]);
    
    %Add randomness to first capillary segment
    if eq(t_init,0)
        length_value = rand(stream)*length_value;
    end
    
    t_end = min(t_init + length_value/v,T);
    length_value = min(length_value,v*(t_end-t_init));
    
    %Compute direction and phase
    delta_pos = length_value * direction;
    pos = pos + delta_pos;
    phase = phase + compute_phase(t_init/T,t_end/T,dt,gG,direction);
    t_init = t_end;
    
    %Store values to stack
    capillary_directions=[capillary_directions;length_value*direction];
end
end

% function [delta_pos] = compute_pos(length,direction)
% direction = delta_pos/norm(delta_pos);
% end

function [phi] = compute_phase(s_1,s_2,dt,gG,v_dir)

%Compute cumulative integral
cum_integral = cumtrapz(dt,gG);

%Compute interval indices
dt_interval = linspace(s_1,s_2,10000);
interval_cum_integral = interp1(dt,cum_integral,dt_interval);

%Compute Integral
integral = trapz(dt_interval,interval_cum_integral);

%Sum phase for each direction
phi=zeros(1,3);

if integral~=0
    phi(1,1) = (integral)*v_dir(1);
    phi(1,2) = (integral)*v_dir(2);
    phi(1,3) = (integral)*v_dir(3);
end
end
