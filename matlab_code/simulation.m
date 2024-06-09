
% Definição dos parâmetros
f = 1; % Frequência da onda (THz)
Delta_x = 10; % FWHM da intensidade transversal (mm)
I0 = 1; % Intensidade central da onda (TW/m^2)
eta_x = 10; % Razão da região de amostragem para FWHM da onda
Nx = 1024; % Número de pontos de amostragem

% Criar uma instância da classe MonoBeam
beam = MonoBeam(f, Delta_x, I0, eta_x, Nx);

% Criar uma máscara de fase (exemplo: máscara de espiral com OAM 1)
M = @(x) exp(1i * angle(x)); % Máscara de fase para OAM 1

% Aplicar a máscara à onda
beam.mask(M);

% Propagar a onda por uma certa distância (exemplo: 10 mm)
beam.propagate(10);

% Plotar a intensidade da onda após a propagação
beam.plotIntensity('Intensidade após a propagação');



