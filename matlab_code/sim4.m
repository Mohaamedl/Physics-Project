
% Parâmetros da grade
N = 1024; % número de pontos na grade
L = 10e-2; % tamanho da grade em metros
dx = L/N; % espaçamento da grade
x = -L/2:dx:L/2-dx; % coordenadas x
y = x; % coordenadas y
[X, Y] = meshgrid(x, y);

% Campo Gaussiano
r = sqrt(X.^2 + Y.^2);
U0 = exp(-(r/w0).^2);

% Máscara de fase espiral (OAM de 1)
l = 2; % ordem do OAM
phi = atan2(Y, X); % ângulo azimutal
spiral_phase_mask = exp(1i * l * phi);
Uin = U0*spiral_phase_mask;
z = 10;
t = ones(1,2,3);
Uout = ang_spec_multi_prop(Uin,w0,0.01,0.1,z,t);
I = abs(Uout).^2;
figure;
imagesc(x, y, I);
colormap('hot');
colorbar;
title('Intensidade do Feixe Gaussiano Após Máscara de Fase Espiral');
xlabel('x (m)');
ylabel('y (m)');
axis square;
%%
% example_pt_source_turb_prop.m

l0 = 0; % inner scale [m]
L0 = inf; % outer scale [m]

zt = [0 z]; % propagation plane locations
Delta_z = zt(2:n) - zt(1:n-1); % propagation distances
% grid spacings
alpha = zt / zt(n);
delta = (1-alpha) * delta1 + alpha * deltan;

% initialize array for phase screens
phz = zeros(N, N, n);
nreals = 20; % number of random realizations
% initialize arrays for propagated fields,
% aperture mask, and MCF
Uout = zeros(N);
mask = circ(xn/D2, yn/D2, 1);
MCF2 = zeros(N);
sg = repmat(sg, [1 1 n]);
for idxreal = 1 : nreals % loop over realizations
    idxreal
    % loop over screens
    for idxscr = 1 : 1 : n
        [phz_lo phz_hi] ...
        = ft_sh_phase_screen ...
        (r0scrn(idxscr), N, delta(idxscr), L0, l0);
        phz(:,:,idxscr) = phz_lo + phz_hi;
    end
    % simulate turbulent propagation
    [xn yn Uout] = ang_spec_multi_prop(pt, wvl, ....
    delta1, deltan, z, sg.*exp(1i*phz));
    % collimate the beam
    Uout = Uout .* exp(-1i*pi/(wvl*R)*(xn.^2+yn.^2));
    % accumulate realizations of the MCF
    %MCF2 = MCF2 + corr2_ft(Uout, Uout, mask, deltan);
end
% modulus of the complex degree of coherence
%MCDOC2 = abs(MCF2) / (MCF2(N/2+1,N/2+1));


function [xn yn Uout] = ang_spec_multi_prop(Uin, wvl, delta1, deltan, z, t)
    % function [xn yn Uout] = ang_spec_multi_prop
    % (Uin, wvl, delta1, deltan, z, t)
    
    N = size(Uin, 1); % number of grid points
    [nx ny] = meshgrid((-N/2 : 1 : N/2 - 1));
    k = 2*pi/wvl; % optical wavevector
    % super-Gaussian absorbing boundary
     nsq = nx.^2 + ny.^2;
     w = 0.47*N;
     sg = exp(-nsq.^8/w^16); clear('nsq', 'w');
    
     z = [0 z]; % propagation plane locations
     n = length(z);
     % propagation distances
     Delta_z = z(2:n) - z(1:n-1);
     % grid spacings
     alpha = z / z(n);
     delta = (1-alpha) * delta1 + alpha * deltan;
     m = delta(2:n) ./ delta(1:n-1);
     x1 = nx * delta(1);
     y1 = ny * delta(1);
    r1sq = x1.^2 + y1.^2;
    Q1 = exp(1i*k/2*(1-m(1))/Delta_z(1)*r1sq);
    Uin = Uin .* Q1 .* t(:,:,1);
    for idx = 1 : n-1
        % spatial frequencies (of i^th plane)
        deltaf = 1 / (N*delta(idx));
        fX = nx * deltaf;
        fY = ny * deltaf;
        fsq = fX.^2 + fY.^2;
        Z = Delta_z(idx); % propagation distance
        % quadratic phase factor
        Q2 = exp(-1i*pi^2*2*Z/m(idx)/k*fsq);
        % compute the propagated field
        Uin = sg .* t(:,:,idx+1) ...
        .* ift2(Q2 ...
        .* ft2(Uin / m(idx), delta(idx)), deltaf);
    end
    % observation-plane coordinates
    xn = nx * delta(n);
    yn = ny * delta(n);
    rnsq = xn.^2 + yn.^2;
    Q3 = exp(1i*k/2*(m(n-1)-1)/(m(n-1)*Z)*rnsq);
    Uout = Q3 .* Uin;
end