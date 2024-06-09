
function mask = azim(l, size, range)
    % Inicializa a máscara como uma matriz de zeros
    mask = zeros(size(1), size(2));
    
    % Define o centro da máscara
    centerX = size(1) / 2;
    centerY = size(2) / 2;
    
    % Define o ângulo de cada fatia da "pizza"
    sliceAngle = 2 * pi / 2^l;
    
    for i = 1:size(1)
        for j = 1:size(2)
            % Calcula o ângulo do ponto atual em relação ao centro
            angle = atan2(i - centerY, j - centerX);
            
            % Ajusta o ângulo para estar entre 0 e 2*pi
            if angle < 0
                angle = angle + 2 * pi;
            end
            
            % Determina a qual fatia da "pizza" o ponto atual pertence
            slice = floor(angle / sliceAngle) + 1;
            
            % Define o valor da célula na máscara
            mask(i, j) = range(1) + mod(slice, 2) * (range(2) - range(1));
        end
    end
end



