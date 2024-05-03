function pattern = hermite_gauss_phase_mask(sz, xmode, ymode, varargin)
    % Generate Hermite-Gaussian mode phase mask
    % Usage:
    %   mask= hermite_gauss_phase_mask([300,300], amode, rmode,...);
    %
    % Parameters:
    %   - sz -- size of the pattern
    %   - xmode -- x mode order
    %   - ymode -- y mode orders
    %
    % Options:
    %   - 'offset'      [x, y] -- offset after applying transformations (default [0,0])
    %   - 'angle'       angle  -- Rotation angle about axis (default: 0 radians)
    %   - 'centre'      [x, y] -- centre location for lens (default: min(sz)/2)
    %   - 'aspect'      aspect -- aspect ratio of lens (default: 1.0)
    %   - 'scale'       scale  -- set the scale (default:min(sz)/10)
    %   - 'range'       [min_val, max_val] -- Range of phase values

    % Set parameters
    p = inputParser;
    setParameters(p, sz);
    lambda = 1555e-9; r = 1e-3;
    s = sqrt(sz(1) * sz(2)) / (2 * max(xmode, ymode));
    %s = sqrt(sz(1) * sz(2)) / max(sz);
    %s = (pi * r^2) / (lambda * (2 * max(xmode, ymode) + 1));
    s = sqrt(sz(1)^2 + sz(2)^2) / (8 * max(xmode, ymode));

    p.addParameter('scale', s);

    p.addParameter('range', [0, 2*pi]); % Default range
    p.parse(varargin{:});

    % Obtain parameters
    c = p.Results.centre;
    ofs = p.Results.offset;
    asp = p.Results.aspect;
    ang = p.Results.angle;
    scale = p.Results.scale;
    range = p.Results.range;
    
    % Generate coordinates
    [xx, yy] = grid2D(sz, 'centre', c, 'offset', ofs, ...
        'angle', ang, 'aspect', asp);

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
    
    pattern = (pattern >= 0);
    pattern = normalize(pattern,'range',range);
end
