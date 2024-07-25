clear all;
clc;

% Carregando as imagens (exemplo)
img1 = imread('1f.png');
img2 = imread('2f.png');
img3 = imread('3f.bmp');
img4 = imread('4f.bmp');
img5 = imread('5f.bmp');

% Agrupando as imagens em uma célula
imgs = {img1, img2, img3, img4, img5};

% Número de imagens
num_imgs = length(imgs);

% Preparando o número de divisões
orders = [1, 2, 3, 4, 5];

% Plotando as imagens em uma linha
figure('Name', 'Beams with Multiple Orders of OAM', 'NumberTitle', 'off');

for i = 1:num_imgs
    % Criando subplot
    subplot(1, num_imgs, i);
    
    % Mostrando a imagem correspondente
    imshow(imgs{i});
    
    % Configurando título acima de cada imagem com o número de divisões
    title(sprintf('Order: %d', orders(i)), 'FontSize', 10);
end

% Ajustando o título principal
sgtitle('Beams with Multiple Orders of OAM', 'FontSize', 14, 'FontWeight', 'bold');
