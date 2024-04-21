function pattern = hermite_gauss_phase_mask(sz, xmode, ymode, varargin)
    % Generate Hermite-Gaussian mode phase mask
    % Usage:
    %   mask= hermite_gauss_phase_mask([300,300], amode, rmode,...);
    %
    % Parameters:
    %   - sz -- size of the pattern
    %   - amode -- azimuthal mode order
    %   - rmode -- radial mode order
    %
    % Options:
    %   - 'centre'      [x, y] -- centre location for lens (default: min(sz)/2)
    %   - 'offset'      [x, y] -- offset after applying transformations
    %   - 'aspect'      aspect -- aspect ratio of lens (default: 1.0)
    %   - 'angle'       angle  -- Rotation angle about axis (radians)
    %   - 'angle_deg'   angle  -- Rotation angle about axis (degrees)
    %   - 'gpuArray'    bool   -- If the result should be a gpuArray
    %   - 'scale'       scale  -- set the scale (default:min(sz)/10)
    

    % Check if the modes xmode and ymode are non-negative integers
    assert(xmode >= 0, 'xmode must be >= 0');
    assert(ymode >= 0, 'ymode must be >= 0');
    assert(floor(xmode) == xmode, 'xmode must be an integer');
    assert(floor(ymode) == ymode, 'ymode must be an integer');
   
    % Parse optional parameters
    p = inputParser;
    setDefaultGridParameters(p, sz);
    %p.addParameter('scale', sqrt(sz(1)^2 +sz(2)^2)/(2*max(xmode, ymode)));
    p.addParameter('scale', min(sz)/10);
    p.parse(varargin{:});

    % Obtain parser parameters
    centre = p.Results.centre;
    offset = p.Results.offset;
    aspect = p.Results.aspect;
    angle = p.Results.angle;
    scale = p.Results.scale;
    % Generate coordinates
    [xx, yy] = grid2D(sz, 'centre', centre, 'offset', offset, ...
        'angle', angle, 'aspect', aspect);

    % Apply scale to coordinates
    xx = xx ./ scale;
    yy = yy ./ scale;

    % Calculate pattern
    xr = linspace(min(xx(:)), max(xx(:)), ceil(sqrt(sz(2)^2 + sz(1)^2)));
    yr = linspace(min(yy(:)), max(yy(:)), ceil(sqrt(sz(2)^2 + sz(1)^2)));

    hx = hermiteH(xmode, cast(xr, 'like', 1));
    hy = hermiteH(ymode, cast(yr, 'like', 1));

    pattern = -ones(sz, 'like', xx);

    for ii = 2:length(xr)
        idx = xx >= xr(ii-1) & xx < xr(ii);
        pattern(idx) = pattern(idx) .* hx(ii-1);
    end
    idx = xx >= xr(end);
    pattern(idx) = pattern(idx) .* hx(end);

    for ii = 2:length(yr)
        idx = yy >= yr(ii-1) & yy < yr(ii);
        pattern(idx) = pattern(idx) .* hy(ii-1);
    end
    idx = yy >= yr(end);
    pattern(idx) = pattern(idx) .* hy(end);

    % Generate phase pattern
    pattern = (pattern >= 0.0) * pi;
end




