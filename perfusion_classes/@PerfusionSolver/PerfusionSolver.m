classdef PerfusionSolver < handle
    
    %TODO: Move some inputs to PerfusionArguments()
    %[k_anisotropy,average_length,sigma_length]
    
    properties(Access = public)
        length_func;
        direction_func;
    end
    
    properties(Access = private)
        inputs=PerfusionArguments();
        average_length=60;
        sigma_length=40;
    end
    
    methods(Access=public)
        function obj = PerfusionSolver(n_particles,varargin)
            obj.inputs.n_particles = n_particles;
            obj.inputs.init_pos = zeros(n_particles,3);
            nInputs = numel(varargin);
            if nInputs > 0
                %Retrieve seed
                obj.inputs.seed = varargin{1};
            else
                fprintf('seed value=1 and n_cores=%d (default).\n',maxNumCompThreads);
            end
        end
        
        function [] = setInitPos(obj,vector)
            obj.inputs.init_pos = vector;
        end
        
        function [] = setGradientProfile(obj,g_profile,T)
            validateattributes(g_profile, {'struct'},{});
            obj.inputs.gradient = g_profile;
            obj.inputs.T = T;
        end
        
        function [real_phase,b_value,Gmax] = computePhase(obj,normalized_phase,varargin)
            gamma = 267.5221900;
            if strcmp(varargin{1},"b-value")
                b_value = varargin{2};
                Gmax = sqrt(b_value/(obj.inputs.T^3*obj.inputs.gradient.a^2*gamma^2));
            elseif strcmp(varargin{1},"Gmax")
                Gmax = varargin{2};
                b_value = gamma^2*Gmax^2*obj.inputs.T^3*obj.inputs.gradient.a^2;
            else
                error("Use either 'b-value' or 'Gmax' a second parameter followed by its value");
            end
            real_phase = obj.inputs.velocity_value * sqrt(b_value*obj.inputs.T) .* normalized_phase;
            
        end
        
        function []=setVelocityValue(obj,velocity_value)
            obj.inputs.velocity_value = velocity_value;
        end
    end
    
    methods(Access=public)
        [phase,capillaries] = computePerfusion(obj,execution_mode);
        [obj] = setCapillaryDistributions(obj,L_distribution,Z_distribution,varargin)
    end
end