% Parâmetros do feixe Gaussiano
lambda = 1555e-9; % Comprimento de onda do laser em metros
w0 = 2e-3; % Raio da cintura do feixe Gaussiano em metros
z = 1; % Distância de propagação em metros

% Parâmetros da máscara de fase Laguerre-Gaussiana
l = 5; % Número quântico azimutal
p = 5; % Número quântico radial


% Criação da grade espacial
N = 1024; % Número de pontos na grade
L = 0.1; % Tamanho da grade em metros
x = linspace(-L/2, L/2, N);
y = linspace(-L/2, L/2, N);
[X, Y] = meshgrid(x, y);
r = sqrt(X.^2 + Y.^2);
phi = atan2(Y, X);
rd = sqrt(N^2 + N^2) / (2 * max(l, p));
% Máscara de fase Laguerre-Gaussiana
Phase = laguerre_gauss_phase_mask([N,N],l,p,'radius',30,'range',[-pi pi]);
%Phase = hermite_gauss_phase_mask([N,N],l,p,'scale',7,'range',[-pi pi]);
% Campo elétrico do feixe Gaussiano
E0 = exp(-(r.^2)/w0^2);

% Aplicação da máscara de fase ao feixe Gaussiano
E_masked = E0 .*exp(1i.*Phase);

% Cálculo da intensidade do feixe resultante após a propagação
E_propagated = fftshift(fft2(fftshift(E_masked)));
I_resultante = abs(E_propagated).^2;

% Visualização da intensidade do feixe resultante
figure;
subplot(1,3,1);
imagesc(x, y, abs(E0).^2);
colormap('hot');
colorbar;
title('Feixe Incidente');
xlabel('x (m)');
ylabel('y (m)');

subplot(1,3,2);
imagesc(x, y, (Phase));
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


% Função Laguerre
function L = laguerreL(p, l, x)
    L = exp(-x/2) .* x.^l .* polyval(LaguerrePoly(p, l), x);
end

% Polinômios de Lagurre associados
function coeffs = LaguerrePoly(p, l)
    coeffs = zeros(1, p+1);
    for k = 0:p
        coeffs(k+1) = ((-1)^k * nchoosek(p+l, p-k)) / factorial(k);
    end
end
