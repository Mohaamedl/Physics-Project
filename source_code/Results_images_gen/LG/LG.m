clear all
clc
lg = {imread("lg_01_100.bmp"), imread("lg_11_80.bmp"), imread("lg_21_85.bmp"), imread("lg_22_110.bmp")};
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

% Campo elétrico do feixe Gaussiano
E0 = exp(-(r.^2)/w0^2);

% Parâmetros para propagação
k = 2 * pi / lambda; % Número de onda
delta = L / N; % Resolução espacial
dz = 0.1;
fX = (-N/2 : N/2-1) / (N*delta);
[x2, y2] = meshgrid(lambda * dz * fX);
h = exp(1i*k/(2*dz)*(x2.^2 + y2.^2)) / (1i*lambda*dz);

% Definindo os modos
modes = [0, 1; 1, 1; 2, 1; 2, 2];
num_modes = size(modes, 1);

% Plotando gráficos para cada tipo de feixe
plot_graphs('Laguerre-Gaussian Beams', @laguerre_gauss_phase_mask, modes, E0, h, delta, x, y, N, lg);

% Função para plotar os gráficos
function plot_graphs(title_text, Phase_func, modes, E0, h, delta, x, y, N, expri)
    figure('Name', title_text, 'NumberTitle', 'off');
    
    num_modes = size(modes, 1);
    for idx = 1:num_modes
        mode1 = modes(idx, 1);
        mode2 = modes(idx, 2);

        % Cálculo da máscara de fase
        Phase = Phase_func([N, N], mode1, mode2, 'radius', 57, "range", [0 2*pi]);

        % Aplicação da máscara de fase
        E_masked = E0 .* exp(1i * Phase);

        % Propagação do feixe
        E_propagated = propagation(E_masked, h, delta);

        % Intensidade simulada
        I_simulated = abs(E_propagated).^2;

        % Plotando os gráficos
        subplot(3, num_modes, idx);
        imagesc(x, y, (Phase));
        axis image;
        
        xlabel('x (m)');
        ylabel('y (m)');
        
        subplot(3, num_modes, idx + num_modes);
        imagesc(x, y, I_simulated);
        axis image;
        
        xlabel('x (m)');
        ylabel('y (m)');

        subplot(3, num_modes, idx + 2 * num_modes);
        resized_img = imresize(expri{idx}, [size(I_simulated, 1) size(I_simulated, 2)]); % Redimensiona a imagem experimental
        imagesc(resized_img); % Utiliza imagesc para manter a escala
        axis image;
        colormap('jet'); % Use gray colormap for experimental images
        
        xlabel('x (m)');
        ylabel('y (m)');
    end
    
    % Adicionando títulos às colunas
    for idx = 1:num_modes
        mode1 = modes(idx, 1);
        mode2 = modes(idx, 2);
        subplot(3, num_modes, idx);
        title(sprintf('Order - (%d,%d)', mode1, mode2));
    end
    
    % Adicionando títulos às linhas
    for row = 1:3
        subplot(3, num_modes, (row-1)*num_modes + 1);
        if row == 1
            ylabel('Phase');
        elseif row == 2
            ylabel('Simulation');
        else
            ylabel('Experimental');
        end
    end
    
    % Ajustando o espaçamento entre os subplots
    ha = get(gcf, 'Children');
    for i = 1:length(ha)
        pos = get(ha(i), 'Position');
        set(ha(i), 'Position', [pos(1), pos(2), pos(3)*0.95, pos(4)*0.95]);
    end
    
    % Ajustando o título principal
    sgtitle(title_text, 'FontSize', 14, 'FontWeight', 'bold', 'Interpreter', 'none');
    
end

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
