img = imbinarize(rgb2gray(im2double(imread('oam_11.bmp'))));
minSize = 35;
img = bwareaopen(img,minSize);
rgs = regionprops(img);
[~,idx] = max([rgs.Area]);
rgs(idx) = [];
center = vertcat(rgs.Centroid);
[c, R] = minboundcircle(center(:,1), center(:,2));
R = round(R);
xc = c(1);yc = c(2);
figure(1)
imshow(img);
hold on
plot(xc, yc, 'b+', 'LineWidth', 2)
images.roi.Circle(gca,'Center',[xc yc],'Radius',R,'Color','g');

%% find the center and radius manually
% figure()
% imshow(img)
% set(gcf, 'pointer', 'arrow');
% zoom on;
% pause() 
% zoom off; 
% [col,row] = ginput(2); % get 2 points, the center and a reference point to get de radius
% zoom out; 
% close()
% c = [col(1) row(1)];
% xc = col(1); yc = row(1);
% 
% 
% 
% 
% 
% R = round(sqrt(abs(col(2)-col(1)).^2 + abs(row(2)-row(1)).^2));
% imshow(img);
% hold on
% plot(xc, yc, 'b+', 'LineWidth', 2)
% images.roi.Circle(gca,'Center',[xc yc],'Radius',R,'Color','g');
% 

%% crop the image


% add padding
img = (rgb2gray(im2double(imread('oam_11.bmp'))));
%a = 30;b = 30; c = s(1); d = s(2);img = imcrop(img,[a b c d]);
s=(size(img));
side = 2*R;
d1 = s(1);
d2 = s(2);
if s(1) <2*R
    d1 = 2*R;
elseif s(2) < 2*R
    d2 = 2*R;
end
img = centralized_mask(img,s(1)+2*R,s(2)+2*R);
figure(2)
imshow(img)

hold on
plot(xc, yc, 'r+', 'LineWidth', 4)

            
%Crop 
img=imcrop(img,[xc yc side side]); 
minSize = 30;
img = imbinarize(img);
img = bwareaopen(img,minSize);
figure(3)
imshow(img)

%%
%preparing for GS
s = max(size(img));
img = ((img - min(img(:))) /(max(img(:)) - min(img(:)))) * 2*pi;
amp = aperture(ones(s,s),4*side/5,s,s);
figure(4)
imshow(amp)
    

%% GS
% 
mask = spiral_mask(s,s,1);
A = mask;
%A = zeros(s,s);
figure(5)
imagesc(A,[0 2*pi])



source= im2double( img);
% 
% A = fftshift(ifft2(fftshift(source)));
iterations = 10;
for It=1:1:iterations
    
    % img plane
    B = amp.*exp(1i.*A);
    BFT = fftshift(fft2(fftshift(B)));
    % s = fftshift(fft2(source));
    % FT plane
    A = angle(BFT)+pi;
    
    C = source.*exp(1i.*A);
    C2 = ifftshift(ifft2(ifftshift(C)));
    % img plane
    A = angle(C2)+pi;


    % B = abs(ring) .* exp(1i*angle(A));
    % C = fftshift(fft2(fftshift(B)));
    % D = abs(source) .* exp(1i*angle(C));
    % A = fftshift(ifft2(fftshift(D)));
 

end

D =angle(exp(1i*A-1i*mask));


figure(6)
imagesc(D)






%%




function [c_mask] = centralized_mask(mask, target_N, target_M)
    [N, M] = size(mask);

    % Calculate the central coordinates for placement
    center_y = floor(target_N / 2);
    center_x = floor(target_M / 2);

    % Calculate the region where the input mask will be placed
    start_y = max(1, center_y - floor(N / 2));
    end_y = min(target_N, start_y + N - 1);
    start_x = max(1, center_x - floor(M / 2));
    end_x = min(target_M, start_x + M - 1);

    % Create a new matrix with zeros
    centralized_phase = zeros(target_N, target_M);

    % Place the input mask at the calculated region
    centralized_phase(start_y:end_y, start_x:end_x) = mask;

    % Resize the matrix to match the target dimensions
    c_mask = centralized_phase(1:target_N, 1:target_M);

    % Convert to binary mask if the input mask was binary
    if islogical(mask)
        c_mask = c_mask > 0.5;
    end
end








function [I] = aperture(I,aperture,w,h)
            
             %Circular aperture
            centerY=round(h/2);
            centerX=round(w/2);
            
            [col, row] = meshgrid(1:w, 1:h);
            circlePixels = (row - centerY).^2 + (col - centerX).^2<= (aperture/2).^2;
            mapa=ones(h,w);
            map=circlePixels.*mapa;
            
            %Aperture addition
            I=map.*I;   
            
end






% credits :
% John D'Errico (2024). A suite of minimal bounding objects
% (https://www.mathworks.com/matlabcentral/fileexchange/34767-a-suite-of-minimal-bounding-objects),
% MATLAB Central File Exchange. Retrieved February 11, 2024.


function [center,radius] = minboundcircle(x,y,hullflag)
% minboundcircle: Compute the minimum radius enclosing circle of a set of (x,y) pairs
% usage: [center,radius] = minboundcircle(x,y,hullflag)
%
% arguments: (input)
%  x,y - vectors of points, describing points in the plane as
%        (x,y) pairs. x and y must be the same size. If x and y
%        are arrays, they will be unrolled to vectors.
%
%  hullflag - boolean flag - allows the user to disable the
%        call to convhulln. This will allow older releases of
%        matlab to use this code, with a possible time penalty.
%        It also allows minboundellipse to call this code
%        efficiently.
% 
%        hullflag = false --> do not use the convex hull
%        hullflag = true  --> use the convex hull for speed
%
%        default: true
%
%
% arguments: (output)
%  center - 1x2 vector, contains the (x,y) coordinates of the
%        center of the minimum radius enclosing circle
%
%  radius - scalar - denotes the radius of the minimum
%        enclosing circle
%
%
% Example usage:
%   x = randn(50000,1);
%   y = randn(50000,1);
%   tic,[c,r] = minboundcircle(x,y);toc
%
%   Elapsed time is 0.171178 seconds.
%
%   c: [-0.2223 0.070526]
%   r: 4.6358
%
%
% See also: minboundrect
%
%
% Author: John D'Errico
% E-mail: woodchips@rochester.rr.com
% Release: 1.0
% Release date: 1/10/07
% default for hullflag
if (nargin<3) || isempty(hullflag)
  hullflag = true;
elseif ~islogical(hullflag) && ~ismember(hullflag,[0 1])
  error 'hullflag must be true or false if provided'
end
% preprocess data
x=x(:);
y=y(:);
% not many error checks to worry about
n = length(x);
if n~=length(y)
  error 'x and y must be the same sizes'
end
% start out with the convex hull of the points to
% reduce the problem dramatically. Note that any
% points in the interior of the convex hull are
% never needed.
if hullflag && (n>3)
  edges = convhulln([x,y]);
  % list of the unique points on the convex hull itself
  % convhulln returns them as edges
  edges = unique(edges(:));
  % exclude those points inside the hull as not relevant
  x = x(edges);
  y = y(edges);
    
end
% now we must find the enclosing circle of those that
% remain.
n = length(x);
% special case small numbers of points. If we trip any
% of these cases, then we are done, so return.
switch n
  case 0
    % empty begets empty
    center = [];
    radius = [];
    return
  case 1
    % with one point, the center has radius zero
    center = [x,y];
    radius = 0;
    return
  case 2
    % only two points. center is at the midpoint
    center = [mean(x),mean(y)];
    radius = norm([x(1),y(1)] - center);
    return
  case 3
    % exactly 3 points
    [center,radius] = enc3(x,y);
    return
end
% more than 3 points.
% Use an active set strategy.
aset = 1:3; % arbitrary, but quite adequate
iset = 4:n;
% pick a tolerance
tol = 10*eps*(max(abs(mean(x) - x)) + max(abs(mean(y) - y)));
% Keep a list of old sets as tried to trap any cycles. we don't need to
% retain a huge list of sets, but only a few of the best ones. Any cycle
% must hit one of these sets. Really, I could have used a smaller list,
% but this is a small enough size that who cares? Almost always we will
% never even fill up this list anyway.
old.sets = NaN(10,3);
old.rads = inf(10,1);
old.centers = NaN(10,2);
flag = true;
while flag
  % have we seen this set before? If so, then we have entered a cycle
  aset = sort(aset);
  if ismember(aset,old.sets,'rows')
    % we have seen it before, so trap out
    center = old.centers(1,:);
    radius = old.radius(1);
    
    % just reset flag then continue, and the while loop will terminate
    flag = false;
    continue
  end
  
  % get the enclosing circle for the current set
  [center,radius] = enc3(x(aset),y(aset));
  
  % is this better than something from the retained sets?
  if radius < old.rads(end)
    old.sets(end,:) = sort(aset);
    old.rads(end) = radius;
    old.centers(end,:) = center;
        
    % sort them in increasing order of the circle radii
    [old.rads,tags] = sort(old.rads,'ascend');
    old.sets = old.sets(tags,:);
    old.centers = old.centers(tags,:);
  end
  
  % are all the inactive set points inside the circle?
  r = sqrt((x(iset) - center(1)).^2 + (y(iset) - center(2)).^2);
  [rmax,k] = max(r);
  if (rmax - radius) <= tol
    % the active set enclosing circle also enclosed
    % all of the inactive points.
    flag = false;
  else
    % it must be true that we can replace one member of aset
    % with iset(k). Which one?
    s1 = [aset([2 3]),iset(k)];
    [c1,r1] = enc3(x(s1),y(s1));
    if (norm(c1 - [x(aset(1)),y(aset(1))]) <= r1)
      center = c1;
      radius = r1;
      
      % update the active/inactive sets
      swap = aset(1);
      aset = [iset(k),aset([2 3])];
      iset(k) = swap;
      
      % bounce out to the while loop
      continue
    end
    s1 = [aset([1 3]),iset(k)];
    [c1,r1] = enc3(x(s1),y(s1));
    if (norm(c1 - [x(aset(2)),y(aset(2))]) <= r1)
      center = c1;
      radius = r1;
      
      % update the active/inactive sets
      swap = aset(2);
      aset = [iset(k),aset([1 3])];
      iset(k) = swap;
      
      % bounce out to the while loop
      continue
    end
    s1 = [aset([1 2]),iset(k)];
    [c1,r1] = enc3(x(s1),y(s1));
    if (norm(c1 - [x(aset(3)),y(aset(3))]) <= r1)
      center = c1;
      radius = r1;
      
      % update the active/inactive sets
      swap = aset(3);
      aset = [iset(k),aset([1 2])];
      iset(k) = swap;
      
      % bounce out to the while loop
      continue
    end
    
    % if we get through to this point, then something went wrong.
    % Active set problem. Increase tol, then try again.
    tol = 2*tol;
    
  end
  
end
end
% =======================================
%  begin subfunction
% =======================================
function [center,radius] = enc3(X,Y)
% minimum radius enclosing circle for exactly 3 points
%
% x, y are 3x1 vectors
% convert to complex
xy = X + sqrt(-1)*Y;
% just in case the points are collinear or nearly so, get
% the interpoint distances, and test the farthest pair
% to see if they work.
Dij = @(XY,i,j) abs(XY(i) - XY(j));
D12 = Dij(xy,1,2);
D13 = Dij(xy,1,3);
D23 = Dij(xy,2,3);
% Find the most distant pair. Test if their circumcircle
% also encloses the third point.
if (D12>=D13) && (D12>=D23)
  center = (xy(1) + xy(2))/2;
  radius = D12/2;
  if abs(center - xy(3)) <= radius
    center = [real(center),imag(center)];
    return
  end
elseif (D13>=D12) && (D13>=D23)
  center = (xy(1) + xy(3))/2;
  radius = D13/2;
  if abs(center - xy(2)) <= radius
    center = [real(center),imag(center)];
    return
  end
elseif (D23>=D12) && (D23>=D13)
  center = (xy(2) + xy(3))/2;
  radius = D23/2;
  if abs(center - xy(1)) <= radius
    center = [real(center),imag(center)];
    return
  end
end
% if we drop down to here, then the points cannot
% be collinear, so the resulting 2x2 linear system
% of equations will not be singular.
A = 2*[X(2)-X(1), Y(2)-Y(1); X(3)-X(1), Y(3)-Y(1)];
rhs = [X(2)^2 - X(1)^2 + Y(2)^2 - Y(1)^2; ...
       X(3)^2 - X(1)^2 + Y(3)^2 - Y(1)^2];
     
center = (A\rhs)';
radius = norm(center - [X(1),Y(1)]);

end




%% mask

function [mask] = spiral_mask(N,M,m)
     x = -N/2:1:(N/2-1); 
    y = -M/2:1:(M/2-1);
    
    % Creates matrices with coordination system
    [x , y] = meshgrid(x,y);
     phi = angle(x+1i*y);               %azimuthal angle
     
     mask=mod(m*phi,2*pi);
    
     
end