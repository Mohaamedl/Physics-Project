clear all
clc

% Carregando as imagens (exemplo)
img1 = imread('4_1024_L1.bmp');
img2 = imread('6_POS.bmp');
img3 = imread('10_POS_L1.bmp');
img4 = imread('14_POS_L1.png');
img5 = imread('16_POS_L1.bmp');
img6 = imread('20_POS.bmp');
img7 = imread('25_POS.bmp');
img8 = imread('30_POS.bmp');

% Agrupando as imagens em uma célula
imgs = {img1, img2, img3, img4, img5, img6, img7, img8};

% Número de imagens
num_imgs = length(imgs);

% Preparando o número de divisões
num_divisoes = [4, 6, 10, 14, 16, 20, 25, 30];

% Plotando as imagens em duas linhas
figure('Name', 'Progress of OAM beam with increasing grating divisions', 'NumberTitle', 'off');

for i = 1:num_imgs
    % Criando subplot
    subplot(2, num_imgs/2, i);
    
    % Mostrando a imagem correspondente
    imshow(imgs{i});
    
    % Configurando título acima de cada imagem com o número de divisões
    title(sprintf('Grating: %d', num_divisoes(i)), 'FontSize', 10);
    
    % Definindo manualmente a posição dos subplots para reduzir o espaço entre as linhas
    if i > num_imgs/2
        pos = get(gca, 'Position');
        set(gca, 'Position', [pos(1), pos(2)-0.07, pos(3), pos(4)]);
    end
end

% Ajustando o título principal
sgtitle('Progress of OAM beam with increasing grating divisions', 'FontSize', 14, 'FontWeight', 'bold');
