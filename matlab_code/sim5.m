% Parâmetros do feixe Gaussiano
lambda = 1550e-9; % Comprimento de onda do laser em metros
w0 = 10e-3; % Raio da cintura do feixe Gaussiano em metros
z = 1; % Distância de propagação em metros

% Parâmetros da máscara de fase Laguerre-Gaussiana
l = 1; % Número quântico azimutal
p = 1; % Número quântico radial

% Criação da grade espacial
N = 2*1024; % Número de pontos na grade
L = 0.2; % Tamanho da grade em metros
x = linspace(-L/2, L/2, N);
y = linspace(-L/2, L/2, N);
[X, Y] = meshgrid(x, y);
r = sqrt(X.^2 + Y.^2);
phi = atan2(Y, X);
rd = sqrt(N^2 + N^2) / (4*(max(l,p)));
% Máscara de fase Laguerre-Gaussiana
Phase = bessel_gauss_phase_mask([N,N], l, 'scale', 7);
Phase = laguerre_gauss_phase_mask([N,N],l,p,'radius',59);
%Phase = (hermite_gauss_phase_mask([N,N],l,p,'scale',59,'range',[0 pi]));

% Campo elétrico do feixe Gaussiano
E0 = exp(-(r.^2)/w0^2);

% Aplicação da máscara de fase ao feixe Gaussiano
E_masked = E0.*exp(1i*Phase);

% Propagação do feixe por convolução
k = 2 * pi / lambda; % Número de onda
delta = L / N; % Resolução espacial

% Fator de fase de propagação 
%H = (exp(1i * k * z *sqrt(1 - (lambda * (r.^2)).^2 / (lambda^2 * z^2))));
fX = (-N/2 : N/2-1) / (L);
% observation-plane coordinates
[x2 y2] = meshgrid(lambda * dz * fX);
clear('fX');
h = exp(1i*k/(2*dz)*(x2.^2+y2.^2)) / (1i*lambda*dz);
% Convolução 
%E_propagated = conv2(E0, E_masked, delta); % Funciona para Hermite

E_propagated = h.*ft2( E_masked, delta); % Funciona para Laguerre e para Bessel
[x2 y2 E_propagated] = one_step_prop(E_masked,lambda,delta,1);
% Cálculo da intensidade do feixe resultante após a propagação
I_resultante = abs(E_propagated).^2;

% Visualização da intensidade do feixe resultante
figure(1);
subplot(1,3,1);
imagesc(x, y, abs(E0).^2);
colormap("hsv");
colorbar;
title('Feixe Incidente');
xlabel('x (m)');
ylabel('y (m)');

subplot(1,3,2);
imagesc(x, y, (Phase));
colormap('hot');
colorbar;
title('Máscara Aplicada');
xlabel('x (m)');
ylabel('y (m)');

subplot(1,3,3);

imagesc(x2, y2, I_resultante);
colormap('parula');
colorbar;
title('Feixe Resultante');
xlabel('x (m)');
ylabel('y (m)');



% Função para convolução 2D
function C = conv2(A, B, delta)
    % function C = myconv2(A, B, delta)
    N = size(A, 1);
    C = ift2(ft2(A, delta) .* ft2(B, delta), 1/(N*delta));
end

function G = ft2(g, delta)
    % function G = ft2(g, delta)
    G = fftshift(fft2(fftshift(g))) * delta^2;
end

function g = ift2(G, delta_f)
    % function g = ift2(G, delta_f)
    N = size(G, 1);
    g = ifftshift(ifft2(ifftshift(G))) * (N * delta_f)^2;
end


function [x2, y2, Uout] = one_step_prop(Uin, wvl, d1, Dz)
% function [x2 y2 Uout] ...
% = one_step_prop(Uin, wvl, d1, Dz)

N = size(Uin, 1); % assume square grid
k = 2*pi/wvl; % optical wavevector
% source-plane coordinates
[x1, y1] = meshgrid((-N/2 : 1 : N/2 - 1) * d1);
 % observation-plane coordinates
 [x2, y2] = meshgrid((-N/2 : N/2-1) / (N*d1)*wvl*Dz);
 % evaluate the Fresnel-Kirchhoff integral
 Uout = 1 / (1i*wvl*Dz) ...
 .* exp(1i * k/(2*Dz) * (x2.^2 + y2.^2)) ...
 .* ft2(Uin .* exp(1i * k/(2*Dz) ...
 * (x1.^2 + y1.^2)), d1);

end