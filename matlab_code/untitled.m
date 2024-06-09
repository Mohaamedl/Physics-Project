





% Carregando a imagem
img = imread("oam_11.bmp");

% Definindo parâmetros
aperture_diameter = 256; % Diâmetro da abertura
padsize = 512; % Tamanho do preenchimento
iterations = 10; % Número de iterações
padding = 512; % Preenchimento adicional

% Chamando a função GS_Correction
correction = GS_Correction(img, aperture_diameter, padsize, iterations, padding);

% Exibindo a correção
imagesc(correction);
colorbar







function correction = GS_Correction(source, aperture_diameter, padsize, iterations, padding)

    % Get position from user input
    [row, col] = get_position();

    % Define half_side
    half_side = round(min(size(source))/4);

    % Crop image
    source = im_crop(source, col(1), row(1), half_side);
    figure(3)
    imshow(source)
    % Convertendo para escala de cinza se for uma imagem colorida (RGB)
    if size(source, 3) == 3
        source = rgb2gray(source);
    end

    % Zero padding of experimental source image
    source = padarray(double(source), [round((padding-size(source,1))/2), round((padding-size(source,2))/2)]);

    % Perfect vortex
    [X1, Y1] = meshgrid(-padsize/2:padsize/2-1, -padsize/2:padsize/2-1);
    phi = angle(X1 + 1i*Y1);
    vortex = mod(11 * phi, 2*pi); % Assuming Topological_charge_Spinner.Value = 11

    % Create the aperture
    amp = aperture(ones(padsize, padsize), aperture_diameter, padsize, padsize);

    % G-S algorithm
    H = vortex;
    E = amp .* exp(1i .* H);
    E_FT = fftshift(fft2(E));
    
    for j = 1:iterations
        E = amp .* exp(1i .* H);
        E_FT = fftshift(fft2(E));
        H = angle(E_FT) + pi;
        disp(size(source));disp(size(H))
        E2 = source .* exp(1i .* H);
        E2 = ifft2(ifftshift(E2));
        H = angle(E2) + pi;
    end
    
    % Preparing correction map
    C = mod(H - vortex, 2*pi);
    disp(size(H))
    C = aperture(C, aperture_diameter, padsize, padsize);
        

    
    disp(size(C))
    % Applying correction to source image
    corrected_image = source .* exp(1i .* C);
    
    % Output the corrected image
    correction = C;
end

% Resto do código permanece inalterado





function [row, col] = get_position()
    img = imread("oam_11.bmp");

    if length(size(img)) ~= 2
        img = im2double(rgb2gray(img));
    end

    %Get position
    figure()
    imshow(img)
    set(gcf, 'pointer', 'arrow');
    zoom on;
    pause() % you can zoom with your mouse and when your image is okay, you press any key
    zoom off; % to escape the zoom mode
    [col, row] = ginput(2);
    zoom out; % go to the original size of your image
    close()
end

function [img] = im_crop(img, center_x, center_y, half_side)
    % Crop image
    x = center_x - half_side;
    y = center_y - half_side;
    side = 2 * half_side;
    img = img(y:y+side-1, x:x+side-1);
end


function [map_aperture] = aperture( aperture, w, h)
    % Circular aperture
    centerY = h / 2;
    centerX = w / 2;
    [columnsInImage, rowsInImage] = meshgrid(1:w, 1:h);
    circlePixels = (rowsInImage - centerY).^2 + (columnsInImage - centerX).^2 <= (aperture/2).^2;
    map_aperture = circlePixels;
end
