function pattern = bessel_gauss_phase_mask(sz, mode, varargin)
% Generates the phase patterns for Bessel beams
% Parameters
%   - sz -- size of the pattern ``[rows, cols]``
%   - mode (integer) -- bessel function mode
% Usage
%   pattern = bessel(sz, mode, ...) generates the phase
%
% Optional named parameters:
%   - 'scale'       scale  -- radial scaling factor for pattern
%   - 'centre'      [x, y] -- centre location for lens (default: sz/2)
%   - 'offset'      [x, y] -- offset after applying transformations
%   - 'type'        type   -- is the lens cylindrical or spherical (1d or 2d)
%   - 'aspect'      aspect -- aspect ratio of lens (default: 1.0)
%   - 'angle'       angle  -- Rotation angle about axis (radians)
%   - 'angle_deg'   angle  -- Rotation angle about axis (degrees)
%   - 'gpuArray'    bool   -- If the result should be a gpuArray


assert(floor(mode) == mode, 'mode must be integer');

p = inputParser;
setDefaultGridParameters(p, sz);
p.addParameter('scale', sqrt(sz(1)^2 + sz(2)^2)/100);
p.parse(varargin{:});
centre = p.Results.centre;
offset = p.Results.offset;
aspect = p.Results.aspect;
ang = p.Results.angle;
% Generate coordinates

[xx,yy] = grid2D(sz, 'centre', centre, 'offset', offset, ...
        'angle', ang, 'aspect', aspect );

r = sqrt(xx.^2 + yy.^2);
% Apply scaling to the coordinates
r = r ./ p.Results.scale;
phi = atan2(yy, xx);
% Calculate the amplitude
amplitude = besselj(mode, r);

% Calculate the phase
pattern = angle(amplitude .* exp(-1i*mode*phi+pi));

pattern = mod(pattern+pi,2*pi);
end