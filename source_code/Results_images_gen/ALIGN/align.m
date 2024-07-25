clear all;
clc;

% Carregando as imagens (exemplo)
img1 = imread('worstH.png');
img2 = imread('interH.png');
img3 = imread('inter2H.png');
img4 = imread('bestH.png');
img5 = imread('worstV.png');
img6 = imread('interV.png');
img7 = imread('inter2V.png');
img8 = imread('bestV.png');

% Agrupando as imagens em uma célula
imgs = {img1, img2, img3, img4, img5, img6, img7, img8};

% Número de imagens
num_imgs = length(imgs);

% Rótulos das imagens com os ângulos de polarização correspondentes
labels = {'H', 'H', 'H', 'H', 'V', 'V', 'V', 'V'};

% Plotando as imagens em duas linhas
figure('Name', 'Wavefront Modulation Progress with Polarization', 'NumberTitle', 'off', 'Color', 'white');

for i = 1:num_imgs
    % Criando subplot
    subplot(2, num_imgs/2, i);
    
    % Mostrando a imagem correspondente
    imshow(imgs{i});
    axis on; % Mostra os eixos para ajustar
    
    % Configurando título acima de cada imagem com o rótulo correspondente
    title(labels{i}, 'FontSize', 14);
    
    % Ajustando os eixos
    set(gca, 'XTick', [], 'YTick', []);
    
    % Definindo manualmente a posição dos subplots para reduzir o espaço entre as linhas
    pos = get(gca, 'Position');
    if i > num_imgs/2
        set(gca, 'Position', [pos(1), pos(2)-0.05, pos(3), pos(4)]);
    end
end

% Ajustando o título principal
sgtitle('Progress of Wavefront Modulation with Polarization', 'FontSize', 16, 'FontWeight', 'bold');

% Melhorando a aparência geral
set(gcf, 'Position', [100, 100, 1200, 600]);
