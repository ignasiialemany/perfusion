function [dt,gG,g_strength,T] = create_seq(config_files)
%% Load sequence
for i = 1:numel(config_files)
    file = config_files{i};
    if ~exist(file, 'file')
        error('Config file "%s" not found', file);
    end
end
configs = cellfun(@yaml.ReadYaml, config_files, 'UniformOutput', false);
config = combine_structs(configs{:});

if isfield(config.sequence, 'type')
    seq_type = config.sequence.type;
    if isfield(config.sequence, seq_type)
        seq_data = config.sequence.(seq_type);
        seq_data.gamma = config.sequence.gamma;
    else
        seq_data = struct();
    end
    seq_Nt = config.sequence.N_t;
    seq_dtmax = cell2mat(config.sequence.dt_max); % one or two
    sequence = MRI.ScanSequence.create(seq_Nt, seq_dtmax, seq_type, seq_data);
else
    % read dt & gG from file, construct sequence directly
    seq = config.sequence.data;
    sequence = MRI.ScanSequence(cell2mat(seq.dt), cell2mat(seq.gG));
end

dt = [0,cumsum(sequence.dt)];
gG = [0,sequence.gG];

g_strength = max(gG);
T = max(dt);

dt = dt./T;
gG = gG./g_strength;
end


function [s_out] = combine_structs(varargin)

% find fieldnames
structures = varargin;
N_structures = numel(structures);
struct_fieldnames = cellfun(@(s) fieldnames(s)', structures, 'UniformOutput', false, 'ErrorHandler', @(err, varargin) '-');
aint_struct = cellfun(@(s) isequal(s, '-'), struct_fieldnames);
if any(aint_struct)
    assert(all(aint_struct), 'field is both struct and value');
    s_out = structures{end}; % take the last value
    return
end
fields = unique([struct_fieldnames{:}]);

% initialise empty output struct
s_out = struct();

for i = 1:numel(fields)
    fieldname = fields{i};
    field_values = cell(1, N_structures);
    for j = 1:numel(structures)
        s = structures{j};
        if isfield(s, fieldname)
            field_values{j} = s.(fieldname);
        end
    end
    field_values = field_values(~cellfun(@isempty, field_values));
    s_out.(fieldname) = combine_structs(field_values{:});
end

end

