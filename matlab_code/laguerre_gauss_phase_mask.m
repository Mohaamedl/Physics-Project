function phase_mask = laguerre_gauss_phase_mask(sz, amode, rmode, varargin)
    % Description:
    % Generates a phase mask for a Gauss-Laguerre (LG) beam.
    % Option:
    %   'radius': Scaling factor for radial mode rings. Default is min(sz)/10.
    % Usage example:
    % phase_mask = generate_lg_phase_mask(sz, amode, rmode);
    %
    % Check mode numbers
    assert(rmode >= 0, 'The radial mode order must be >= 0');
    assert(floor(rmode) == rmode, 'The radial mode order must be an integer');
    assert(floor(amode) == amode, 'The azimuthal mode order must be an integer');

    % Parse the input parameters
    p = inputParser;
    setDefaultGridParameters(p,sz);  
    p.addParameter('radius', min(sz)/10);
    p.addParameter('p0', 1.0);
    p.parse(varargin{:});

    % Obtain parser parameters
    c = p.Results.centre;
    ofs = p.Results.offset;
    asp = p.Results.aspect;
    ang = p.Results.angle;
    
    % Generate coordinates
    [xx, yy] = grid2D(sz, 'centre', c, 'offset', ofs, ...
        'aspect', asp, 'angle',ang);

    % Calculate the azimuthal part of the pattern
    phase = amode .* atan2(yy, xx);

    % Calculate the Laguerre polynomials in the radial direction
    rho = sqrt(xx.^2 + yy.^2);
    Lpoly_rho = fastLaguerre(rmode, abs(amode), (rho./p.Results.radius).^2);

    % Calculate the radial part of the phase
    phase = phase + (sign(Lpoly_rho)>0)*pi;

    phase_mask = mod(phase, 2*pi);
    
end
