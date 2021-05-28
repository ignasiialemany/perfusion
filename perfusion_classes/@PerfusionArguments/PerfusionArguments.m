classdef PerfusionArguments
    
    properties(Access = public)
        n_particles = 1000;
        init_pos;
        velocity_value = 0.4;
        gradient = struct();
        T=1;
        G=1;
        seed=1;
    end
    
    methods
        
        %Number of Particles
        function obj = set.n_particles(obj,n_particles)
            validateattributes( n_particles, { 'numeric' }, { 'scalar', 'positive' } );
            obj.n_particles = n_particles;
        end
        
        %Initial Positions
        function obj = set.init_pos(obj,init_pos)
            validateattributes( init_pos, { 'numeric' }, { '2d' } );
            obj.init_pos = init_pos;
        end
        
        %Blood flow velocity value
        function obj = set.velocity_value(obj,velocity_value)
            validatattributes( velocity_value, { 'numeric' }, {'finite','scalar'});
            obj.velocity_value = velocity_value;
        end
        
        %Load and compute Gradient Profile
        function obj=set.gradient(obj,gradient)
            obj.gradient.gG = gradient.gG/max(gradient.gG);
            obj.gradient.dt = gradient.dt/max(gradient.dt);
            obj.gradient.a = sqrt(trapz(obj.gradient.dt,cumtrapz(obj.gradient.dt,obj.gradient.gG).^2));
        end
    end
end