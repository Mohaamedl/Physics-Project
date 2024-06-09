% Parâmetros do feixe Gaussiano
lambda = 1555e-9; % Comprimento de onda do laser em metros
w0 = 10e-3; % Raio da cintura do feixe Gaussiano em metros
z = 0.01; % Distância de propagação em metros

% Parâmetros da máscara de fase Laguerre-Gaussiana
l = 2; % Número quântico azimutal
p = 2; % Número quântico radial

% Criação da grade espacial
N = 1024; % Número de pontos na grade
L = 0.1; % Tamanho da grade em metros
x = linspace(-L/2, L/2, N);
y = linspace(-L/2, L/2, N);
[X, Y] = meshgrid(x, y);
r = sqrt(X.^2 + Y.^2);
phi = atan2(Y, X);

% Máscara de fase Laguerre-Gaussiana
Phase = laguerre_gauss_phase_mask([N,N],l,p);

% Campo elétrico do feixe Gaussiano
E0 = exp(-(r.^2)/w0^2);

% Aplicação da máscara de fase ao feixe Gaussiano
E_masked = E0 .* exp(1i * Phase);

% Propagação do feixe por convolução
k = 2 * pi / lambda; % Número de onda
delta = L / N; % Resolução espacial

% Fator de fase de propagação (transformada de Fourier da função de propagação)
H = exp(1i * k * z * sqrt(1 - (lambda * (X.^2 + Y.^2)).^2 / (lambda^2 * z^2)));

% Convolução no domínio da frequência
E_propagated = myconv2(E_masked, H, delta); 

% Cálculo da intensidade do feixe resultante após a propagação
I_resultante = abs(E_propagated).^2;

% Visualização da intensidade do feixe resultante
figure(1);
subplot(1,3,1);
imagesc(x, y, abs(E0).^2);
colormap('hot');
colorbar;
title('Feixe Incidente');
xlabel('x (m)');
ylabel('y (m)');

subplot(1,3,2);
imagesc(x, y, angle(Phase));
colormap('hot');
colorbar;
title('Máscara Aplicada');
xlabel('x (m)');
ylabel('y (m)');

subplot(1,3,3);
imagesc(x, y, I_resultante);
colormap('hot');
colorbar;
title('Feixe Resultante');
xlabel('x (m)');
ylabel('y (m)');

% Função que cria a máscara de fase Laguerre-Gaussiana
function Phase = laguerre_gauss_phase_mask(size, l, p)
    % Parâmetros da máscara de fase
    radius = max([l, p]) * 2; % Ajuste automático do raio
    range = [-pi pi]; % Faixa padrão de fase
    
    % Criação da grade espacial
    x = linspace(-radius, radius, size(1));
    y = linspace(-radius, radius, size(2));
    [X, Y] = meshgrid(x, y);
    [theta, rho] = cart
end