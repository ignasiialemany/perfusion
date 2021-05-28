function [tensor, RMS] = process_signal(data, sequence)

% combine IC and EC data from simulation
phase = [data.phase_ECS; data.phase_ICS];
displ = [data.displacement_ECS; data.displacement_ICS];

% calculate diffusion tensor
t = cumsum(sequence.dt);
k = cumtrapz(t, sequence.gG);
bvalue = trapz(t,k.^2);

tensor = process_phase(phase, bvalue);
tensor = MRI.DiffusionTensor(tensor);

% calculate bulk diffusivity
T = sum(sequence.dt); % NOT equal to Delta, because we don't have that displ data
RMS = process_displacement(displ, T);

end

function [tensor] = process_phase(phase, bvalue)
% get diffusion tensor

% gradient sampling directions
directions = [1,  1,  0;
              1, -1,  0;
              1,  0,  1;
              1,  0, -1;
              0,  1,  1;
              0,  1, -1];

ndirs = size(directions, 1);
b_matrix = zeros(ndirs, 6);
signal_ratio = zeros(ndirs, 1);
for i = 1:ndirs
    dir = directions(i, :); % [Gx, Gy, Gz]
    % b-matrix (b_xx, b_yy, b_zz, b_xy, b_xz, b_yz)
    b_matrix(i, :) = [dir([1, 2, 3]).^2, 2*dir([1, 1, 2]).*dir([2, 3, 3])] * bvalue;
    % attenuation vector
    phi = sum(dir .* phase, 2); % combine components
    signal_ratio(i) = abs(mean(exp(-1i*phi), 1));
end

% least squares solution
A = b_matrix;
b = log(signal_ratio);
x = -lscov(A, b); % solve

% assign tensor
tensor = zeros(3, 3);
tensor(tril(true(3))) = x([1, 4, 5, 2, 6, 3]); % assign the right spaces
tensor(triu(true(3), +1)) = tensor(tril(true(3), -1)); % mirror over diagonal

end

function [RMS] = process_displacement(displacement, T)
% get diffusivity from RMS displacement

dim = size(displacement, 2);
RMS.MD = mean(sum(displacement.^2, 2))/(2*dim*T); % Cartesian distance
RMS.Dx = mean(displacement(:, 1).^2)/(2*T); % distance in x only
RMS.Dy = mean(displacement(:, 2).^2)/(2*T); % distance in y only
RMS.Dz = mean(displacement(:, 3).^2)/(2*T); % distance in z only

end
