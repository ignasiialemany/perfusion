classdef ScanSequence

    methods
        function obj = ScanSequence(dt, gG)
            N = numel(dt);
            obj.dt(1, 1:N) = dt;
            obj.gG(:, 1:N) = gG;
        end
    end

    properties(SetAccess=private)
        dt(1,:)
        gG(1,:) %TODO: support 1D & 3D gradient vector
    end
    methods
        function dt_n = get_dt(obj, n)
            dt_n = obj.dt(1, n);
        end
        function gG_n = get_gG(obj, n)
            gG_n = obj.gG(:, n);
        end
    end

    properties(Dependent)
        NT
        bvalue
    end
    methods
        function [NT] = get.NT(obj)
            NT = length(obj.dt);
        end
        function [b] = get.bvalue(obj)
            t = cumsum(obj.dt);
            k = cumtrapz(t, obj.gG);
            b = trapz(t, k.^2);
        end
    end

    methods(Static)
        [sequence] = create(NT, dt_max, SeqName, varargin)
    end

end
