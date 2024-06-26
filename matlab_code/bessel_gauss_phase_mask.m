function pattern = bessel_gauss_phase_mask(sz, mode, varargin)
    % Generates the phase patterns for Bessel-Gaussian beams
    % Parameters
    %   - sz -- size of the pattern ``[rows, cols]``
    %   - mode (integer) -- bessel function mode
    % Usage
    %   pattern = bessel_gauss_phase_mask(sz, mode, ...);
    %
    % Optional named parameters:
    %   - 'scale'       scale  -- radial scaling factor for pattern
    %   - 'centre'      [x, y] -- centre location for lens (default: sz/2)
    %   - 'offset'      [x, y] -- offset after applying transformations
    %   - 'aspect'      aspect -- aspect ratio of lens (default: 1.0)
    %   - 'angle'       angle  -- Rotation angle about axis (radians)
    
    % Set parameters
    p = inputParser;
    setParameters(p, sz);
    p.addParameter('scale', sqrt(sz(1)^2 + sz(2)^2)/100);
    p.addParameter('range', [0 2*pi]);
    p.parse(varargin{:});


    c = p.Results.centre;
    ofs = p.Results.offset;
    asp = p.Results.aspect;
    ang = p.Results.angle;
    range = p.Results.range;
    % Generate coordinates
    [xx,yy] = grid2D(sz, 'centre', c, 'offset', ofs, ...
            'angle', ang, 'aspect', asp );
    
    r = sqrt(xx.^2 + yy.^2);
    % Apply scaling to the coordinates
    r = r ./ p.Results.scale;
    phi = atan2(yy, xx);
    % Calculate the amplitude
    amplitude = besselj(mode, r);
    
    % Calculate the phase
    pattern = angle(amplitude .* exp(-1i*mode*phi+pi)) + pi;
    phase = mod(pattern+pi,2*pi);
    phase_range = range(2) - range(1);

    % Normalize phase to the specified range
    phase_diff = range(1) - min(phase(:));
    phase_mask = phase + phase_diff;
    
    % Ensure phase is within the specified range
    pattern = mod(phase_mask - range(1), phase_range) + range(1);
end