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

% Campo elétrico do feixe Gaussiano
E0 = exp(-(r.^2)/w0^2);

% Parâmetros para propagação
k = 2 * pi / lambda; % Número de onda
delta = L / N; % Resolução espacial
dz = 0.1;
fX = (-N/2 : N/2-1) / (N*delta);
[x2, y2] = meshgrid(lambda * dz * fX);
h = exp(1i*k/(2*dz)*(x2.^2 + y2.^2)) / (1i*lambda*dz);
clear('fX');


% Figura 1: Laguerre-Gauss
figure;
for l_laguerre = 0:3
    for p_laguerre = 0:3
        % Cálculo das máscaras de fase para Laguerre-Gauss
        Phase_laguerre = laguerre_gauss_phase_mask([N, N], l_laguerre, p_laguerre, 'radius', 27);
        % Aplicação das máscaras de fase aos feixes Gaussianos
        E_masked_laguerre = E0.*exp(1i*Phase_laguerre);
        % Propagação do feixe
        E_propagated_laguerre = propagation(E_masked_laguerre, h.*E_masked_laguerre, delta);
        % Cálculo das intensidades do feixe resultante
        I_resultante_laguerre = abs(E_propagated_laguerre).^2;
        % Visualização do resultado
        subplot(4, 4, l_laguerre*4 + p_laguerre + 1);
        imagesc(x, y, I_resultante_laguerre);
        colormap("jet");
        colorbar;
        title(sprintf("LG_{(%d,%d)}", l_laguerre, p_laguerre));
        xlabel('x (m)');
        ylabel('y (m)');
    end
end
sgtitle('Intensidade dos Feixes Propagados - Laguerre-Gauss');

% Figura 2: Hermite-Gauss
figure;
for x_hermite = 0:3
    for y_hermite = 0:3
        % Cálculo das máscaras de fase para Hermite-Gauss
        Phase_hermite = hermite_gauss_phase_mask([N, N], x_hermite, y_hermite, 'scale', 26, 'range', [0 pi]);
        % Aplicação das máscaras de fase aos feixes Gaussianos
        E_masked_hermite = E0.*exp(1i*Phase_hermite);
        % Propagação do feixe
        E_propagated_hermite = propagation(E_masked_hermite, h.*E_masked_hermite, delta);
        % Cálculo das intensidades do feixe resultante
        I_resultante_hermite = abs(E_propagated_hermite).^2;
        % Visualização do resultado
        subplot(4, 4, x_hermite*4 + y_hermite + 1);
        imagesc(x, y, I_resultante_hermite);
        colormap("jet");
        colorbar;
        title(sprintf("HG_{(%d,%d)}", x_hermite, y_hermite));
        xlabel('x (m)');
        ylabel('y (m)');
    end
end
sgtitle('Intensidade dos Feixes Propagados - Hermite-Gauss');

% Figura 3: Bessel-Gauss
figure;
for m_bessel = -2:2
    % Cálculo das máscaras de fase para Bessel-Gauss
    Phase_bessel = bessel_gauss_phase_mask([N, N], m_bessel, 'scale', 7);
    % Aplicação das máscaras de fase aos feixes Gaussianos
    E_masked_bessel = E0.*exp(1i*Phase_bessel);
    % Propagação do feixe
    E_propagated_bessel = propagation(E_masked_bessel, h.*E_masked_bessel, delta);
    % Cálculo das intensidades do feixe resultante
    I_resultante_bessel = abs(E_propagated_bessel).^2;
    % Visualização do resultado
    subplot(1, 5, m_bessel + 3);
    imagesc(x, y, I_resultante_bessel);
    colormap("jet");
    colorbar;
    title(sprintf("Bessel, m = %d", m_bessel));
    xlabel('x (m)');
    ylabel('y (m)');
end
sgtitle('Intensidade dos Feixes Propagados - Bessel-Gauss');
% Função para convolução 2D
function C = conv2D(A, B, delta)
    N = size(A, 1);
    C = ift2(ft2(A, delta) .* ft2(B, delta), 1/(N*delta));
end

function G = ft2(g, delta)
    G = fftshift(fft2(fftshift(g))) * delta^2;
end

function g = ift2(G, delta_f)
    N = size(G, 1);
    g = ifftshift(ifft2(ifftshift(G))) * (N * delta_f)^2;
end

% Função para propagação
function E_propagated = propagation(E_masked, h, delta)
    E_propagated = h .* conv2D(E_masked, E_masked, delta);
end
