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
    s = sqrt(sz(1)^2 + sz(2)^2) / (8 * max(xmode, ymode));

    p.addParameter('scale', min(sz)/10);
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
    hx = hermite(xmode, cast(xr, 'like', 1));
    hy = hermite(ymode, cast(yr, 'like', 1));
    
    % Interpolate Hermite polynomials
    hx_interp = interp1(xr, hx, xx, 'nearest', 'extrap');
    hy_interp = interp1(yr, hy, yy, 'nearest', 'extrap');
    
    % Generate pattern using vectorized operations
    pattern = prod(cat(3, hx_interp, hy_interp), 3);
    
    % Generate phase pattern
    pattern = ~(pattern >= 0);
    pattern = normalize(pattern, 'range', range);
end


% Avan Suinesiaputra (2024). Hermite polynomials (https://www.mathworks.com/matlabcentral/fileexchange/27746-hermite-polynomials),
% MATLAB Central File Exchange. Retrieved May 13, 2024.
function h = hermite(n,x)
% HERMITE: compute the Hermite polynomials.
% 
%   h = hermite(n)
%   h = hermite(n,x)
% 
% Inputs:
%   - n is the order of the Hermite polynomial (n>=0).
%   - x is (optional) values to be evaluated on the resulting Hermite
%     polynomial function.
% 
% There are two possible outputs:
% 1. If x is omitted then h is an array with (n+1) elements that contains
%    coefficients of each Hermite polynomial term.
%    E.g. calling h = hermite(3)
%    will result h = [8 0 -12 0], i.e. 8x^3 - 12x
% 
% 2. If x is given, then h = Hn(x) and h is the same size of x.
%    E.g., H2(x) = 4x^2 - 2
%    calling h = hermite(2,[0 1 2])
%    will result h = [-2 2 14]
% 
% More information:
% - about the Hermite polynomial: http://mathworld.wolfram.com/HermitePolynomial.html
% - some examples of this function:
% http://suinotes.wordpress.com/2010/05/26/hermite-polynomials-with-matlab/
% 
% Authors: Avan Suinesiaputra (avan.sp@gmail.com)
%          Fadillah Z Tala    (fadil.tala@gmail.com)
% rev.
% 26/05/2010 - first creation.
%            - bug fixed: error when hermite(0,x) is called (x isn't empty)
% 24/09/2010 - bug fixed: the size of x does match with y in line 50.
%              (thanks to Shiguo Peng)
% check n
if( n<0 ), error('The order of Hermite polynomial must be greater than or equal to 0.'); end
% again check n is an integer
if( 0~=n-fix(n) ), error('The order of Hermite polynomial must be an integer.'); end
% call the hermite recursive function.
h = hermite_rec(n);
% evaluate the hermite polynomial function, given x
if( nargin==2 )
    y = h(end) * ones(size(x));
    p = 1;
    for i=length(h)-1:-1:1
        y = y + h(i) * x.^p;
        p = p+1;
    end
    
    % restore the shape of y, the same as x
    h = reshape(y,size(x));
end
function h = hermite_rec(n)
% This is the reccurence construction of a Hermite polynomial, i.e.:
%   H0(x) = 1
%   H1(x) = 2x
%   H[n+1](x) = 2x Hn(x) - 2n H[n-1](x)
if( 0==n ), h = 1;
elseif( 1==n ), h = [2 0];
else
    
    h1 = zeros(1,n+1);
    h1(1:n) = 2*hermite_rec(n-1);
    
    h2 = zeros(1,n+1);
    h2(3:end) = 2*(n-1)*hermite_rec(n-2);
    
    h = h1 - h2;
    
end
end
end