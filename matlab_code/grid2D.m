function [xx, yy] = grid2D(sz, varargin)
    % Generate a grid of points with some options
    % Options:
    %   - 'centre'      [x, y] -- centre location for lens (default: min(sz)/2)
    %   - 'offset'      [x, y] -- offset after applying transformations
    %   - 'aspect'      aspect -- aspect ratio of lens (default: 1.0)
    %   - 'angle'       angle  -- Rotation angle about axis (radians)
    %   - 'angle_deg'   angle  -- Rotation angle about axis (degrees)
    %   - 'gpuArray'    bool   -- If the result should be a gpuArray
    

    p = inputParser;
    setDefaultGridParameters(p, sz);
    parse(p, varargin{:});

    if p.Results.gpuArray
        [xx, yy] = meshgrid(gpuArray(1:sz(2)), gpuArray(1:sz(1)));
    else
        [xx, yy] = meshgrid(1:sz(2), 1:sz(1));
    end

    xx = xx - p.Results.centre(1);
    yy = yy - p.Results.centre(2);

    
    angle = p.Results.angle;
    if isempty(angle)
        angle = deg2rad(p.Results.angle_deg);
    end
    xx = cos(angle) .* xx - sin(angle) .* yy;
    yy = sin(angle) .* xx + cos(angle) .* yy;
    yy = yy * p.Results.aspect;

    xx = xx - p.Results.offset(1);
    yy = yy - p.Results.offset(2);
end