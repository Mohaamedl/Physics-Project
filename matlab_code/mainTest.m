clear all
close all
clc
%%
% Example of use LG

figure(1)

for i = 0:4
    for ii = 0:4
        amode = i;
        rmode = ii;
        mask = laguerre_gauss_phase_mask([256, 256], amode, rmode,'aspect',1, ...
           'offset',[0,0],'angle',0,'radius',25);
        subplot(5, 5, i * 5 + ii + 1)
        imagesc(mask, [0 2*pi])
        title(sprintf("LG_{(%d,%d)}", rmode, amode));
        colormap hot
        colorbar
    end
end


%%

% Example Usage HG:
figure(2)
set(gcf,'Position',[500 200 800 600])
for i = 0:4
    for ii = 0:4
        amode = i;
        rmode = ii;
        mask = hermite_gauss_phase_mask([256, 256], amode, rmode,'aspect',1,'offset',[0,0],'scale',26);
        subplot(5, 5, i * 5 + ii + 1)
        imagesc(mask, [0 pi])
        title(sprintf("HG_{(%d,%d)}", amode, rmode));
        colorbar
    end
end

%%
% Example of use Bessel

% Define the parameters of the Bessel mode
sz = [512, 512];  % Size of the mask

figure(3)
for i = -2:2
    mode = i;
    phase_mask = bessel_gauss_phase_mask(sz, mode, 'scale', 5);
    subplot(2, 3, i + 3)
    imagesc(phase_mask, [0 2*pi]);
    colorbar
    colormap hot
    title(sprintf("BG_{%d}", mode));
end



