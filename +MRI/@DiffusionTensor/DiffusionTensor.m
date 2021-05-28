% bare minimum classdef required to read DiffusionTensor from .mat file
classdef DiffusionTensor
    properties
        tensor(3, 3)
    end
    methods
        function obj = DiffusionTensor(tensor)
            obj.tensor = tensor;
        end
    end
end
