function phase_mask = laguerre_gauss_phase_mask(sz, amode, rmode, varargin)
    % Generates the phase patterns for Laguerre-Gaussian beams
    % Parameters
    %   - sz -- size of the pattern [rows, cols]
    %   - amode -- azimuthal mode order
    %   - rmode -- radial mode orders 
    % Usage
    %   pattern = laguerre_gauss_phase_mask(sz, amode, rmode,...);
    %
    % Optional named parameters:
    %   - 'scale'       scale  -- radial scaling factor for pattern
    %   - 'centre'      [x, y] -- centre location for lens (default: sz/2)
    %   - 'offset'      [x, y] -- offset after applying transformations
    %   - 'aspect'      aspect -- aspect ratio of lens (default: 1.0)
    %   - 'angle'       angle  -- Rotation angle about axis (radians)

    % Set parameters
    p = inputParser;
    setParameters(p,sz);  
    p.addParameter('radius', min(sz)/10);
    p.parse(varargin{:});

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
