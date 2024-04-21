function setParameters(p, sz, varargin)
    % Set default values for grid parameters 

    assert(numel(sz) == 2, 'sz must be a vector of 2 elements');

    parameters = {'offset', 'angle','centre',  'aspect'};

    ip = inputParser;

    ip.addParameter('offset', [0, 0]);
    ip.addParameter('angle', 0);
    ip.addParameter('centre', (sz)/2);
    ip.addParameter('aspect', 1.0);
    ip.addParameter('skip', {});
    ip.parse(varargin{:});

    for pp = 1:length(parameters)
        if ~any(strcmpi(parameters{pp}, ip.Results.skip))
            p.addParameter(parameters{pp}, ip.Results.(parameters{pp}));
        end
    end
end
