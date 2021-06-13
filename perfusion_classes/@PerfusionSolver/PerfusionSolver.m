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
        k_anisotropy=10;
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
        
        function [obj] = setCapillary(obj,length_distribution,zenith_distribution)
            
            
           
            
            
        end
        function [obj] = setLengthDistribution(obj,distribution,varargin)
            if strcmp(distribution,'Weibull')
                obj.length_func = @(x) weibull(x,obj.average_length,obj.sigma_length);
            elseif strcmp(distribution,'Constant')
                N=varargin{1};
                obj.length_func = @(x) (obj.inputs.T * obj.inputs.velocity_value)/N;
            else
                fprintf("This length distribution has not been implemented yet");
            end
        end
        
        function [obj] = setAnglesDistributions(obj,distribution)
            if strcmp(distribution,'Anisotropic')
                obj.zenith_func = @(x,mean) vonmises(x,mean,obj.k_anisotropy);
                obj.azimuth_func = @(x) 2*pi*x;
            elseif strcmp(distribution,'Isotropic')
                obj.zenith_func = @(x,mean) acos(-1+2*x);
                obj.azimuth_func = @(x) 2*pi*x;
            else
                fprintf("This angle distribution has not been implemented yet");
            end
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
    end
    
    methods(Access=public)
        [phase,capillaries] = computePerfusion(obj,execution_mode);
        [obj] = setCapillaryDistributions(obj,L_distribution,Z_distribution)
    end
end