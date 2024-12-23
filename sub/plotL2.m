function plotL2(L2, dataPath, outFile)
fontName = machine_prefs;
fileName = strrep(outFile,'.mat','.png');
plotPath = fullfile(dataPath,'Plots/L2/');
if ~isfolder(plotPath)
    mkdir(plotPath)
end

fh = figure;
% for s =1:length(sensorList)
subplot(4,1,1)
plot(L2.wavelength,L2.Irradiance.Es, 'LineWidth', 2)
hold on
plot(L2.wavelength, L2.Irradiance.Es + L2.Irradiance.Es_sd,'--k')
plot(L2.wavelength, L2.Irradiance.Es - L2.Irradiance.Es_sd,'--k')
text(0.82,0.9,sprintf('Hourly Spec Filter: %.1f%%',...
    100 * L2.QC.specFilter.Es_Hourly_Nremoved/L2.QC.specFilter.Es_Hourly_Nstart),...
    'Units','normalized', 'FontSize', 12)
grid on
title(sprintf('Es: %s', strrep(fileName,'_',' ')))
ylabel('[W/m^2/nm]')
set(gca, 'FontName', fontName, 'FontSize', 18)
xlim([335 935])

subplot(4,1,2)
plot(L2.wavelength,L2.Radiance.Li, 'LineWidth', 2)
hold on
plot(L2.wavelength, L2.Radiance.Li + L2.Radiance.Li_sd,'--k')
plot(L2.wavelength, L2.Radiance.Li - L2.Radiance.Li_sd,'--k')
text(0.82,0.9,sprintf('Hourly Spec Filter: %.1f%%',...
    100 * L2.QC.specFilter.Li_Hourly_Nremoved/L2.QC.specFilter.Li_Hourly_Nstart),...
    'Units','normalized', 'FontSize', 12)
grid on
title('L_i')
ylabel('[W/m^2/sr/nm]')
set(gca, 'FontName', fontName, 'FontSize', 18)
xlim([335 935])

subplot(4,1,3)
plot(L2.wavelength,L2.Radiance.Lt, 'LineWidth', 2)
hold on
plot(L2.wavelength, L2.Radiance.Lt + L2.Radiance.Lt_sd,'--k')
plot(L2.wavelength, L2.Radiance.Lt - L2.Radiance.Lt_sd,'--k')

plot(L2.wavelength_rho,L2.Radiance.Lw, 'LineWidth', 2)
hold on
plot(L2.wavelength_rho, L2.Radiance.Lw + L2.Radiance.Lw_unc,'--k')
plot(L2.wavelength_rho, L2.Radiance.Lw - L2.Radiance.Lw_unc,'--k')

text(0.82,0.9,sprintf('Hourly Spec Filter: %.1f%%',...
    100 * L2.QC.specFilter.Lt_Hourly_Nremoved/L2.QC.specFilter.Lt_Hourly_Nstart),...
    'Units','normalized', 'FontSize', 12)
grid on
title('L_t & L_w')
ylabel('[W/m^2/sr/nm]')
set(gca, 'FontName', fontName, 'FontSize', 18)
xlim([335 935])
plot([335 935],[0 0],'-k', 'LineWidth', 1.5)

subplot(4,1,4)
plot(L2.wavelength_rho,L2.Reflectance.Rrs, 'LineWidth', 2)
hold on
plot(L2.wavelength_rho, L2.Reflectance.Rrs + L2.Reflectance.Rrs_unc,'--k')
plot(L2.wavelength_rho, L2.Reflectance.Rrs - L2.Reflectance.Rrs_unc,'--k')
% text(0.75,0.75,sprintf('AVW: %.1f',L2.Reflectance.AVW(1)),...
%     'fontname',fontName,'fontSize',14,'Units','normalized')
grid on
title('R_{rs}')
ylabel('[1/sr]')
xlim([335 935])
plot([335 935],[0 0],'-k', 'LineWidth', 1.5)
ylim([-2e-4 max(max(L2.Reflectance.Rrs + L2.Reflectance.Rrs_unc))])
text(0.01,0.9,sprintf('Negative Rrs Ens Filtered: %d',...
    L2.QC.negRrsRemoved),...
    'Units','normalized', 'FontSize', 12)
for i=1:size(L2.Reflectance.Rrs,1)
    text(0.72,1.05-0.1*i,sprintf('AVW: %.0f QWIP: %.3f N: %d',...
        L2.Reflectance.AVW(i), L2.Reflectance.QWIP(i), ...
        L2.QC.glitterStats(i).finalN), ...
        'Units','normalized', 'FontSize', 12)
end
set(gca, 'FontName', fontName, 'FontSize', 18)
set(gcf,"Position",[1862        -177        1119        1166])

exportgraphics(fh,sprintf('%s%s',plotPath,fileName))