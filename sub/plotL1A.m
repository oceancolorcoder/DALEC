%% Figures
function plotL1A(sensorList,L1A,fileName,plotPath,flags)
for s = 1:length(sensorList)
    fh = figure;
    subplot(2,1,1)
    text(0.05, 0.95, sprintf('Start: %s End: %s',min(L1A.(sensorList{s}).datetime), max(L1A.(sensorList{s}).datetime)) )
    text(0.05, 0.85, sprintf('Lat: %.3f +/- %.3f',mean(L1A.(sensorList{s}).lat), std(L1A.(sensorList{s}).lat)) )
    text(0.05, 0.75, sprintf('Lon: %.3f +/- %.3f',mean(L1A.(sensorList{s}).lon), std(L1A.(sensorList{s}).lon)) )
    text(0.05, 0.65, sprintf('SZA: %.1f +/- %.1f',mean(L1A.(sensorList{s}).sza), std(L1A.(sensorList{s}).sza)) )
    text(0.05, 0.55, sprintf('RelAz: %.1f +/- %.1f',mean(L1A.(sensorList{s}).relAz), std(L1A.(sensorList{s}).relAz)) )
    text(0.05, 0.45, sprintf('Pitch: %.1f +/- %.1f',mean(L1A.(sensorList{s}).pitch), std(L1A.(sensorList{s}).pitch)) )
    text(0.05, 0.35, sprintf('Roll: %.1f +/- %.1f',mean(L1A.(sensorList{s}).roll), std(L1A.(sensorList{s}).roll)) )
    text(0.05, 0.25, sprintf('Temp: %.1f +/- %.1f',mean(L1A.(sensorList{s}).temp), std(L1A.(sensorList{s}).temp)) )
    % text(0.05, 0.15, sprintf('QFlag (already filtered): %d of %d', sum(L1A.(sensorList{s}).qFlag), length(L1A.(sensorList{s}).qFlag)) )
    text(0.05, 0.15, sprintf('IntTime: %.1f +/- %.1f',mean(L1A.(sensorList{s}).intTime), std(L1A.(sensorList{s}).intTime)) )
    title(strrep(fileName,'_',' '))
    axis off

    subplot(2,1,2)
    ph1 = plot( mean(L1A.(sensorList{s}).light),'.k');
    hold on
    ph2 = plot( mean(L1A.(sensorList{s}).light) + std(L1A.(sensorList{s}).light),'--k');
    plot(mean(L1A.(sensorList{s}).light) - std(L1A.(sensorList{s}).light),'--k')
    ph3 = plot( mean(repmat(L1A.(sensorList{s}).dark, 1,size(L1A.(sensorList{s}).light, 2))) ,'.r');
    text(0.75,0.9,sprintf('Start N: %d',flags.(sensorList{s}).startN),'units','normalized')
    text(0.75,0.8,sprintf('GPS: %d (%.1f%%)',flags.(sensorList{s}).GPS,100*(flags.(sensorList{s}).GPS/flags.(sensorList{s}).startN)),'units','normalized')
    text(0.75,0.7,sprintf('qFlag: %d (%.1f%%)',flags.(sensorList{s}).qFlag,100*(flags.(sensorList{s}).qFlag/flags.(sensorList{s}).startN)),'units','normalized')
    text(0.75,0.6,sprintf('relAz: %d (%.1f%%)',flags.(sensorList{s}).relAz,100*(flags.(sensorList{s}).relAz/flags.(sensorList{s}).startN)),'units','normalized')
    text(0.75,0.5,sprintf('sza: %d (%.1f%%)',flags.(sensorList{s}).sza,100*(flags.(sensorList{s}).sza/flags.(sensorList{s}).startN)),'units','normalized')
    text(0.75,0.4,sprintf('tilt: %d (%.1f%%)',flags.(sensorList{s}).tilt,100*(flags.(sensorList{s}).tilt/flags.(sensorList{s}).startN)),'units','normalized')
    text(0.75,0.3,sprintf('nans: %d (%.1f%%)',flags.(sensorList{s}).nans,100*(flags.(sensorList{s}).nans/flags.(sensorList{s}).startN)),'units','normalized')
    text(0.75,0.2,sprintf('spec: %d (%.1f%%)',flags.(sensorList{s}).spec,100*(flags.(sensorList{s}).spec/flags.(sensorList{s}).startN)),'units','normalized')
    title(sprintf('Final N: %d',size(L1A.(sensorList{s}).light,1)))

    legend([ph1 ph2 ph3],{sensorList{s},'STD','Darks'}, 'Location','northwest')
    grid on
    xlabel('pixel')
    ylabel('counts [x10^4]')
    title(sensorList{s})

    fprintf('Saving %s%s_L1A_%s.png\n',plotPath,fileName,(sensorList{s}))
    exportgraphics(fh,sprintf('%s%s_L1A_%s.png',plotPath,fileName,(sensorList{s})))
    close
end
