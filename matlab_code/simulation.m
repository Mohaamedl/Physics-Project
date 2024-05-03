% Parâmetros do feixe
lambda = 632.8e-9; % Comprimento de onda em metros
w0 = 1e-3; % Raio do feixe no waist (cintura) em metros
z = 0.1; % Distância de propagação em metros

% Coordenadas transversais
x = linspace(-5e-3, 5e-3, 1000);
y = x;
[X, Y] = meshgrid(x, y);
r = sqrt(X.^2 + Y.^2);

% Cálculo do campo do feixe Gaussiano
k = 2*pi/lambda; % Número de onda
zR = pi*w0^2/lambda; % Raio de curvatura de Rayleigh
wz = w0*sqrt(1+(z/zR)^2); % Raio do feixe em z
Rz = z*(1+(zR/z)^2); % Raio de curvatura da frente de onda em z
phi = atan(z/zR); % Fase de Gouy
E = (w0/wz)*exp(-r.^2/wz^2).*exp(-1i*(k*z + k*r.^2/(2*Rz) - phi));

% Visualização do feixe Gaussiano
figure(1);
imagesc(x, y, abs(E).^2);
title('Intensidade do Feixe Gaussiano');
xlabel('x (m)');
ylabel('y (m)');
colorbar;

% Parâmetros do feixe Laguerre-Gaussiano
l = 3; % Ordem azimutal
p = 2; % Ordem radial

% Cálculo do campo do feixe Laguerre-Gauss
rho = sqrt(X.^2 + Y.^2);
theta = atan2(Y,X);
Lpl = @(p,l,rho) ((rho.^l).*exp(-rho.^2/2).*fastLaguerre(p,l,rho.^2)); % Polinômio de Laguerre generalizado
E_LG = (w0./wz).*rho.^l .* Lpl(p,l,2*rho.^2/wz^2) .* exp(-1i*l*theta) .* exp(-rho.^2./wz.^2) .* exp(-1i*k*z - 1i*k*rho.^2/(2*Rz) + 1i*l*phi);

% Visualização do feixe Laguerre-Gauss
figure(2);
imagesc(x, y, abs(E_LG).^2);
title('Intensidade do Feixe Laguerre-Gauss');
xlabel('x (m)');
ylabel('y (m)');
colorbar;
%%
% Parâmetros do feixe Hermite-Gaussiano
lambda = 632.8e-9; % Comprimento de onda do feixe (em metros)
w0 = 1e-3; % Raio da cintura do feixe no plano z=0 (em metros)
z = 0; % Posição ao longo do eixo z (em metros)
n_ref = 1; % Índice de refração do meio
k = 2*pi*n_ref/lambda; % Número de onda

% Parâmetros dos polinômios de Hermite-Gauss
m = 0; % Ordem do polinômio de Hermite (modo transversal)
n = 0; % Ordem do polinômio de Hermite (modo transversal)

% Grade de coordenadas transversais
x = linspace(-5*w0, 5*w0, 100);
y = linspace(-5*w0, 5*w0, 100);
[X, Y] = meshgrid(x, y);

% Raio do feixe em z
wz = w0 * sqrt(1 + (lambda*z/(pi*w0^2))^2);

% Curvatura da frente de onda em z
Rz = z * (1 + (pi*w0^2/(lambda*z))^2);

% Fase Gouy em z
Xi = atan(z * lambda/(pi*w0^2));

% Campo elétrico do feixe Hermite-Gaussiano
E = ((w0/wz) * exp(-((X.^2+Y.^2)/wz^2)) .* ...
     polyval(hermiteH(m, sqrt(2)*X/wz), sqrt(2)*X/wz).' .* ...
     polyval(hermiteH(n, sqrt(2)*Y/wz), sqrt(2)*Y/wz) .* ...
     exp(-1i*k*z - 1i*k*(X.^2+Y.^2)/(2*Rz) + 1i*(m+n+1)*Xi));

% Intensidade do feixe
I = abs(E).^2;

% Plot da intensidade do feixe
figure(3);
imagesc(x, y, I);
colormap('hot');
colorbar;
xlabel('x (m)');
ylabel('y (m)');
title('Intensidade do Feixe Hermite-Gaussiano');
