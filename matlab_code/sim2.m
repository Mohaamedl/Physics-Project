% Parâmetros do feixe Gaussiano
lambda = 1555.8e-9; % comprimento de onda do laser em metros
w0 = 1e-3; % raio da cintura do feixe Gaussiano em metros
z = 1; % distância de propagação em metros

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

% Aplicação da máscara de fase ao campo Gaussiano
U = U0 .* spiral_phase_mask.*exp(1i*lambda/(2*pi)*z);

% Simulação de propagação
U_prop = fftshift(fft2(fftshift(U)));

% Intensidade resultante
I = abs(U_prop).^2;

% Visualização
figure;
imagesc(x, y, I);
colormap('hot');
colorbar;
title('Intensidade do Feixe Gaussiano Após Máscara de Fase Espiral');
xlabel('x (m)');
ylabel('y (m)');
axis square;
