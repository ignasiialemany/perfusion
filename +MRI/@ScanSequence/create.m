function [sequence] = create(NT, dt_max, SeqName, varargin)
% make a sequence object from the given parameters
%   NT : number of time steps (target)
%   dt_max : maximum time step
%       (1, 2) = (dt_free, dt_grad) - in case of real sequence
%       scalar = dt - in case of dummy sequence
%   SeqName : sequence name {PGSE, MCSE/M2SE, STEAM}
%   varargin : parameters specific to the given sequence with SeqName

switch upper(SeqName)
    case {'PGSE', 'STEAM'} % follow the same rules
        [Gmax, alpha90, alphaRO, epsilon, Delta, delta, gamma] = parse(varargin, ...
            {'Gmax', 'alpha90', 'alphaRO', 'epsilon', 'Delta', 'delta', 'gamma'});
        durations = [alpha90, ...
                     epsilon, delta, epsilon, ...
                     Delta-(2*epsilon+delta), ...
                     epsilon, delta, epsilon, ...
                     alphaRO];
        ids = [0, 1, 2, 3, 0, -1, -2, -3, 0];
        [dt, gG] = discretize(durations, ids, NT, dt_max, Gmax);
    case {'MCSE', 'M2SE'} % these two names are synonyms
        [Gmax, alpha90, alphaRO, epsilon, delta1, delta2, gamma] = parse(varargin, ...
            {'Gmax', 'alpha90', 'alphaRO', 'epsilon', 'delta1', 'delta2', 'gamma'});
        del1 = delta1+2*epsilon;
        del2 = delta2+2*epsilon;
        Delta = (del2*(-2*del1+epsilon) + del1*epsilon)/(del1-del2);
        durations = [alpha90, ...
                     epsilon, delta1, epsilon, epsilon, delta2, epsilon, ...
                     Delta-(del1+del2), ...
                     epsilon, delta1, epsilon, epsilon, delta2, epsilon, ...
                     alphaRO];
        ids = [0, -1, -2, -3, 1, 2, 3, 0, -1, -2, -3, 1, 2, 3, 0];
        [dt, gG] = discretize(durations, ids, NT, dt_max, Gmax);
    otherwise % dummy sequence
        gamma = 1;
        dt = ones(1, NT)*dt_max(1);
        gG = zeros(1, NT);
end

gG = gG*gamma; % apply gamma

sequence = MRI.ScanSequence(dt, gG);

end

function [varargout] = parse(args, names)

if isscalar(args) && isstruct(args{1})
    argstruct = args{1};
    N = length(names);
    varargout = cell(1, N);
    for i = 1:N
        varargout{i} = argstruct.(names{i});
    end
else
    varargout = args;
end

end

function [dt, gG] = discretize(durations, ids, NT, dt_max, Gmax)
% discretize the sequence
%   durations := length of intervals
%   ids := designtaion of interval
%       0 - flat (free, gradient OFF)
%       % gradients: +/- sign of id represents sign of Gmax at the end of the gradient
%       1 - (+/-) gradient ramp-up
%       2 - (+/-) gradient flat
%       3 - (+/-) gradient ramp-down

% calculate the target time step dt
dt_aim = sum(durations)/NT; % to get approximately NT steps
% enforce limits
dt_free = min(dt_aim, dt_max(1));
dt_grad = min(dt_aim, dt_max(2));

% calculate the number of steps per interval, rounded up based on target dt
Nt_intervals(ids==0) = ceil(durations(ids==0)/dt_free);
Nt_intervals(ids~=0) = ceil(durations(ids~=0)/dt_grad);

dt = [];
gG = [];
for i = 1:numel(durations)
    % time steps
    Nt_i = Nt_intervals(i);
    dt_i = durations(i)/Nt_i;
    dt = [dt, repmat(dt_i, [1, Nt_i])];
    % gradient
    id_i = ids(i);
    switch abs(id_i) % use absolute value for switch and use sign(id_i) inside cases
        case 0 % flat, gradient off
            gA = 0;
            gB = gA;
        case 1 % ramp-up gradient
            gA = 0;
            gB = sign(id_i)*Gmax;
        case 2 % flat, gradient on
            gA = sign(id_i)*Gmax;
            gB = gA;
        case 3 % ramp-down gradient
            gA = sign(id_i)*Gmax;
            gB = 0;
    end
    gvals = linspace(gA, gB, Nt_i+1);
    gG = [gG, gvals(1:end-1)];
end

end
