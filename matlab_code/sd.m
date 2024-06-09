clear all
clc



% Parâmetros comuns
lambda = 1555e-9; % Comprimento de onda do laser em metros
w0 = 10e-3; % Raio da cintura do feixe Gaussiano em metros
z = 1; % Distância de propagação em metros
% Criação da grade espacial
N = 2*1024; % Número de pontos na grade
L = 0.2; % Tamanho da grade em metros
x = linspace(-L/2, L/2, N);
y = linspace(-L/2, L/2, N);
[X, Y] = meshgrid(x, y);
r = sqrt(X.^2 + Y.^2);
phi = atan2(Y, X);
rd = sqrt(N^2 + N^2) / (4);

% Parâmetros para Bessel
m_bessel = 3; % Modo


% Parâmetros para Laguerre-Gauss
l_laguerre = 2; % Número quântico azimutal 
p_laguerre = 2; % Número quântico radial

% Parâmetros para Hermite
x_hermite = 2; % Ordem x
y_hermite = 2; % Ordem y

% Cálculo das máscaras de fase para cada caso
Phase_bessel = bessel_gauss_phase_mask([N, N], m_bessel, 'scale', 7);

Phase_laguerre = laguerre_gauss_phase_mask([N, N], l_laguerre, p_laguerre, 'radius', 57);

Phase_hermite = hermite_gauss_phase_mask([N, N], x_hermite, y_hermite, 'scale', 27, 'range', [0 pi]);


% Campo elétrico do feixe Gaussiano
E0 = exp(-(r.^2)/w0^2);

% Aplicação das máscaras de fase aos feixes Gaussianos
E_masked_bessel = E0.*exp(1i*Phase_bessel);
E_masked_laguerre = E0.*exp(1i*Phase_laguerre);
E_masked_hermite = E0.*E0.*exp(1i*Phase_hermite);

% Propagação dos feixes por convolução

k = 2 * pi / lambda; % Número de onda
delta = L / N; % Resolução espacial
% Convolução 
dz = 0.1;
fX = (-N/2 : N/2-1) / (N*delta);
% observation-plane coordinates
[x2, y2] = meshgrid(lambda * dz * fX);
h = exp(1i*k/(2*dz)*(x2.^2+y2.^2)) / (1i*lambda*dz);
clear('fX');
H = exp(1i*k*dz)/(1i*lambda*dz)*exp(1i*k/2 *dz*(x2^2+y2^2));
isequal(H,h)
E_propagated_bessel = h.*conv2D((E_masked_bessel), (E_masked_bessel), delta);
E_propagated_laguerre = h.* conv2D(E_masked_laguerre, (E_masked_laguerre), delta);
E_propagated_hermite = h.*conv2D(E_masked_hermite, abs(E_masked_hermite), delta);

% Cálculo das intensidades dos feixes resultantes
I_resultante_bessel = abs(E_propagated_bessel).^2;
I_resultante_laguerre = abs(E_propagated_laguerre).^2;
I_resultante_hermite = abs(E_propagated_hermite).^2;
clear('h')


% Visualização dos resultados
figure(1);

% Para Bessel
subplot(3, 3, 1);
imagesc(x, y, abs(E0).^2);
colormap("jet");
colorbar;
title('Feixe Incidente');
xlabel('x (m)');
ylabel('y (m)');

subplot(3, 3, 2);
imagesc(x, y, (Phase_bessel));
colorbar;
title(['Máscara Aplicada - Bessel, m = ', num2str(m_bessel)]);
xlabel('x (m)');
ylabel('y (m)');

subplot(3, 3, 3);
imagesc(x, y, I_resultante_bessel);
colorbar;
title('Feixe Resultante');
xlabel('x (m)');
ylabel('y (m)');

% Para Laguerre-Gauss
subplot(3, 3, 4);
imagesc(x, y, abs(E0).^2);
colorbar;
title('Feixe Incidente');
xlabel('x (m)');
ylabel('y (m)');

subplot(3, 3, 5);
imagesc(x, y, (Phase_laguerre));
colorbar;
title(sprintf("Máscara Aplicada - LG_{(%d,%d)}", l_laguerre, p_laguerre));
xlabel('x (m)');
ylabel('y (m)');

subplot(3, 3, 6);
imagesc(x, y, I_resultante_laguerre);
colorbar;
title('Feixe Resultante');
xlabel('x (m)');
ylabel('y (m)');

%Para Hermite
subplot(3, 3, 7);
imagesc(x, y, abs(E0).^2);
colorbar;
title('Feixe Incidente');
xlabel('x (m)');
ylabel('y (m)');

subplot(3, 3, 8);
imagesc(x, y, (Phase_hermite));
colorbar;
title(sprintf("Máscara Aplicada - HG_{(%d,%d)}", x_hermite, y_hermite));

xlabel('x (m)');
ylabel('y (m)');

subplot(3, 3, 9);
imagesc(x, y, I_resultante_hermite);
colorbar;
title('Feixe Resultante');
xlabel('x (m)');
ylabel('y (m)');
% figure(4)
% im = (im2gray(I_resultante_laguerre));
% imagesc(im)
% im = GerchbergSaxton(im);
% imagesc(angle(im))



%%






% Campo elétrico do feixe Gaussiano

x = linspace(-L/2, L/2, N);
y = linspace(-L/2, L/2, N);
[X, Y] = meshgrid(x, y);
r = sqrt(X.^2 + Y.^2);
E0 = exp(-(r.^2)/w0^2);
E_masked= E0.*exp(1i*phase);
k = 2 * pi / lambda; % Número de onda
delta = L / N; % Resolução espacial
fX = (-N/2 : N/2-1) / (N*delta);
% observation-plane coordinates
[x2, y2] = meshgrid(lambda * dz * fX);
h = exp(1i*k/(2*dz)*(x2.^2+y2.^2)) / (1i*lambda*dz);
clear('fX');
e = propagation(Phase_hermite);

function E_propagated = propagation(phase)
    E_masked= E0.*exp(1i*phase);
    E_propagated = h.*conv2D((E_masked), (E_masked), delta);
end



























% Função para convolução 2D
function C = conv2D(A, B, delta)
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







%% 
% trash
% figure(1)
% for i=1:length(z)
%     E_propagated_bessel(:,:,i) = h.*conv2(E_masked_bessel, E_masked_bessel, delta);
%     E_propagated_laguerre(:,:,i) =h.* conv2(E_masked_laguerre, E_masked_laguerre, delta);
%     E_propagated_hermite(:,:,i) = h.*conv2(E_masked_hermite, abs(E_masked_hermite), delta);
%     imagesc(x, y,abs(E_propagated_bessel(:,:,i).^2) );
%     hold on
%     pause(0.1)
% end
% hold off
