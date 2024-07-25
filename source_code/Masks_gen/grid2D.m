function [xx, yy] = grid2D(sz, varargin)
    % Generate a grid of points with some options:

    %   - 'offset'      [x, y] -- offset after applying transformations
    %   (default [0,0])
    %   - 'angle'       angle  -- Rotation angle about axis (default: 0 radians)
    %   - 'centre'      [x, y] -- centre location for lens (default: min(sz)/2)
    %   - 'aspect'      aspect -- aspect ratio of lens (default: 1.0)

    % Set parameteres
    p = inputParser;
    setParameters(p, sz);
    parse(p, varargin{:});
    % Generate coordinates
    [xx, yy] = meshgrid(1:sz(2), 1:sz(1));
    % Centralize
    xx = xx - p.Results.centre(1);
    yy = yy - p.Results.centre(2);

    % Rotate
    ang = p.Results.angle;
    
    xx = cos(ang) .* xx - sin(ang) .* yy;
    yy = sin(ang) .* xx + cos(ang) .* yy;
    % Reshape yy
    yy = yy * p.Results.aspect;
    % Translate
    xx = xx - p.Results.offset(1);
    yy = yy - p.Results.offset(2);
end