
clear all

close all
% Defina o alvo
guesss = laguerre_gauss_phase_mask([512,512],1,0)/(2*pi);

% Defina o número de iterações
num_iterations = 100;

n=512;
x = linspace(-10,10,n);
y = linspace(-10,10,n);
[X,Y] = meshgrid(x,y);
x0 = 0;     		% center
y0 = 0;     		% center
sigma = 2; 			% beam waist
A = 1;      		% power
targ = A  * exp(-((X-x0).^2 + (Y-y0).^2)./(2*sigma^2));
[X,Y] = meshgrid(x,y);
[Phi,R] = cart2pol(X,Y);
R0 = 6; %ring radius
W = 2; %ring width
ring = exp(-((R-R0).^2)./W.^2);
dinput_intensity = double(rgb2gray(imread("oam_11.bmp")));
input = dinput_intensity(110:621,270:781);
minSize = 6000;
input = bwareaopen(imbinarize(input),minSize);
input = ((input - min(input(:))) /(max(input(:)) - min(input(:)))) ;
guess = input;
targ = ring;
figure(1)
imshow(targ)
figure(2)
% Execute o algoritmo Gerchberg-Saxton
mask = gerchberg_saxton(ring,guess, num_iterations);
imshow(angle(mask))


function final_guess = gerchberg_saxton(target, initial_guess, num_iterations)
    % Executa o algoritmo Gerchberg-Saxton para estimar a fase do padrão de holograma.
    %
    % Parâmetros:
    %   - target: padrão alvo complexo
    %   - initial_guess: palpite inicial para o padrão de holograma
    %   - num_iterations: número de iterações do algoritmo
    %
    % Saída:
    %   - final_guess: melhor estimativa da fase do padrão de holograma após
    %     o número especificado de iterações
    
    % Verifique se o palpite inicial foi fornecido
    if nargin < 3
        error('Número insuficiente de argumentos. Por favor, forneça o alvo, o palpite inicial e o número de iterações.');
    end
    
    % Inicialize o palpite atual como o palpite inicial
    current_guess = initial_guess;
    
    % Loop de iterações
    for i = 1:num_iterations
        % Calcula o padrão gerado a partir do palpite atual
        generated_pattern = fftshift(fft2(ifftshift(current_guess)));
        
        % Calcula o padrão alvo de amplitude com a fase do padrão gerado
        target_amplitude = abs(target) .* exp(1i * angle(generated_pattern));
        
        % Calcula o novo palpite usando a transformada inversa
        current_guess = fftshift(ifft2(ifftshift(target_amplitude)));
        
        
    end
    
    % A saída final_guess contém a melhor estimativa da fase do padrão de holograma
    final_guess = current_guess;
end
