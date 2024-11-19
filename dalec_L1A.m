% Read raw dalec file into arrays

% Screen for NaNs in GPS Lat/Lon
% Screen for QFlag "This flag for your current DALEC firmware version (v3.5-96-g130bc2f) relates
%   to the GearPos movement. 0 = stationary and 1=moving during the spectrometer integration. This
%   column will contain extra bits of info with the next firmware update and will let you know of the details when the time comes."
% Screen for relAz
% Screen for SZA
% Screen for tilt
% Screen for spectral outliers and NaNs
%
% Input: Raw text files from DALEC acquisition software (duration unclear)
% Output: L1A matlab structures
%
% D. Aurin NASA/GSFC November 2024

%% Setup
wipe
rawPath = '~/Projects/HyperPACE/field_data/DALEC/VIIRS2024/raw';
L1Apath = '~/Projects/HyperPACE/field_data/DALEC/VIIRS2024/L1A';

fileList = dir(fullfile(rawPath,'*.TXT'));
sensorList = {'Ed','Lsky','Lu'};
minMaxRelAz = [88 137]; % Zhang 2017: 90*, Mobley 1999: 135, Zibordi 2009 (and IOCCG Protocols): 90
minMaxSZA = [15 65];    %e.g. 20: Zhang 2017, depends on wind, 60:Brewin 2016
maxTilt = 5;            % 2-5 deg. IOCCG Draft Protocols
minSpectra = 5;         % If any sensor has fewer than this # of spectra, don't save this file
sigmaLight.Ed = 5;      % For filtering spectral outliers
sigmaLight.Lsky = 8;
sigmaLight.Lu = 3;
% fRange = [350, 900];    % Spectral range for spectral filter
fRange = [5 170];    % Wavelength is not yet defined. Based on cals for pixels.

if ~isfolder(L1Apath)
    mkdir(L1Apath);
end

makePlots = 1;

%% Process L1A
% There are different numbers of records for each instrument

for i=1:length(fileList)
    failFlag = 0;
    fpf = fullfile(fileList(i).folder,fileList(i).name);
    % Detect end of header
    fid = fopen(fpf,"r");
    for l=1:100
        tline = fgetl(fid);
        if contains(tline,"OUTPUT FORMAT")
            break
        end
    end
    fclose(fid);

    opts = detectImportOptions(fpf,...
        "ExpectedNumVariables",212,...
        "NumHeaderLines",l,...
        "Delimiter",",");

    fprintf('Processing: %s\n',fpf)
    T = readtable(fpf,opts);

    % There is a mismatch between the field names and the datasets starting
    % with UTCtimestamp and GPSfix.
    sensor = T.UTCtimestamp;

    for s = 1:length(sensorList)
        Idx = strcmp(sensor,sensorList{s});
        L1A.(sensorList{s}).datetime = datetime(T.GPSfix(Idx),'Format','yyyy-MM-dd''T''HH:mm:ss.SSSZ','TimeZone','UTC');
        L1A.(sensorList{s}).lat = T.Lat(Idx);
        L1A.(sensorList{s}).lon = T.Lon(Idx);
        L1A.(sensorList{s}).sza = T.SolarZenith(Idx);
        L1A.(sensorList{s}).relAz = abs(T.RelAz(Idx)); % Ignore sign for now.
        L1A.(sensorList{s}).pitch = T.Pitch(Idx);
        L1A.(sensorList{s}).roll = T.Roll(Idx);
        L1A.(sensorList{s}).temp = T.DetectorTemp(Idx);
        L1A.(sensorList{s}).qFlag = T.Qflag(Idx);
        L1A.(sensorList{s}).intTime = T.Inttime(Idx);
        L1A.(sensorList{s}).dark = T.DarkCounts(Idx);
        L1A.(sensorList{s}).light = table2array(T(Idx,23:212));

        fieldNames = fieldnames(L1A.(sensorList{s}));

        N = length(L1A.(sensorList{s}).lat);

        % GPS
        badGPS = isnan(L1A.(sensorList{s}).lat) | isnan(L1A.(sensorList{s}).lon);
        fprintf('%s records of %d removed for bad GPS: %d\n', sensorList{s}, length(badGPS), sum(badGPS))
        for f=1:length(fieldNames)
            if contains(fieldNames{f},'light')
                L1A.(sensorList{s}).(fieldNames{f})(badGPS,:) = [];
            else
                L1A.(sensorList{s}).(fieldNames{f})(badGPS) = [];
            end
        end

        % Quality Flag
        badQ = logical(L1A.(sensorList{s}).qFlag);
        fprintf('%s records of %d removed for bad QFlag: %d\n', sensorList{s}, length(badQ), sum(badQ))
        for f=1:length(fieldNames)
            if contains(fieldNames{f},'light')
                L1A.(sensorList{s}).(fieldNames{f})(badQ,:) = [];
            else
                L1A.(sensorList{s}).(fieldNames{f})(badQ) = [];
            end
        end

        % Relative Azimuth
        badRelAz = L1A.(sensorList{s}).relAz < minMaxRelAz(1) | ...
            L1A.(sensorList{s}).relAz > minMaxRelAz(2)| ...
            isnan(L1A.(sensorList{s}).relAz);
        fprintf('%s records of %d removed for bad RelAz: %d\n', sensorList{s}, length(badRelAz), sum(badRelAz))
        for f=1:length(fieldNames)
            if contains(fieldNames{f},'light')
                L1A.(sensorList{s}).(fieldNames{f})(badRelAz,:) = [];
            else
                L1A.(sensorList{s}).(fieldNames{f})(badRelAz) = [];
            end
        end

        % SZA
        badSZA = L1A.(sensorList{s}).sza < minMaxSZA(1) | ...
            L1A.(sensorList{s}).sza > minMaxSZA(2)| ...
            isnan(L1A.(sensorList{s}).sza);
        fprintf('%s records of %d removed for bad SZA: %d\n', sensorList{s}, length(badSZA), sum(badSZA))
        for f=1:length(fieldNames)
            if contains(fieldNames{f},'light')
                L1A.(sensorList{s}).(fieldNames{f})(badSZA,:) = [];
            else
                L1A.(sensorList{s}).(fieldNames{f})(badSZA) = [];
            end
        end

        % Tilt
        badTilt = L1A.(sensorList{s}).pitch > maxTilt | ...
            L1A.(sensorList{s}).roll > maxTilt | ...
            isnan(L1A.(sensorList{s}).roll);
        fprintf('%s records of %d removed for bad tilt: %d\n', sensorList{s}, length(badTilt), sum(badTilt))
        for f=1:length(fieldNames)
            if contains(fieldNames{f},'light')
                L1A.(sensorList{s}).(fieldNames{f})(badTilt,:) = [];
            else
                L1A.(sensorList{s}).(fieldNames{f})(badTilt) = [];
            end
        end

        % NaNs 
        badNans = isnan(L1A.(sensorList{s}).light);
        [badNans,~] = find(badNans);
        badNans = unique(badNans);
        fprintf('%s records of %d removed for NaNs: %d\n', sensorList{s}, length(badNans), length(badNans))
        for f=1:length(fieldNames)
            if contains(fieldNames{f},'light')
                L1A.(sensorList{s}).(fieldNames{f})(badNans,:) = [];
            else
                L1A.(sensorList{s}).(fieldNames{f})(badNans) = [];
            end
        end

        % Spectral Outliers        
        if length(L1A.(sensorList{s}).datetime) >= minSpectra         
            subLight = L1A.(sensorList{s}).light(:,fRange(1):fRange(2));
            peakLight = max(L1A.(sensorList{s}).light(:,fRange(1):fRange(2)),[],2);
            normLight = subLight ./ peakLight;

            ave = mean(normLight);
            stD = std(normLight);
            sigma = sigmaLight.(sensorList{s});
            flag = normLight < (ave - sigma*stD) |...
                normLight > (ave + sigma*stD);
            [badSpec,~] = find(flag);
            badSpec = unique(badSpec);
            fprintf('%s records of %d removed for spectral outliers: %d\n', sensorList{s}, length(flag), length(badSpec))
            for f=1:length(fieldNames)
                if contains(fieldNames{f},'light')
                    L1A.(sensorList{s}).(fieldNames{f})(badSpec,:) = [];
                else
                    L1A.(sensorList{s}).(fieldNames{f})(badSpec) = [];
                end
            end
            clear normLight ave stD flag
        end
        flags.(sensorList{s}).startN = N;
        flags.(sensorList{s}).GPS = sum(badGPS);
        flags.(sensorList{s}).qFlag = sum(badQ);
        flags.(sensorList{s}).relAz = sum(badRelAz);
        flags.(sensorList{s}).sza = sum(badSZA);
        flags.(sensorList{s}).tilt = sum(badTilt);
        flags.(sensorList{s}).nans = length(badNans);
        flags.(sensorList{s}).spec = length(badSpec);


        if length(L1A.(sensorList{s}).datetime) < minSpectra
            disp('**************************************************')
            fprintf('Fewer than %d %s spectra remaining. FAIL. DO NOT SAVE!\n',minSpectra, sensorList{s})
            disp('**************************************************')
            failFlag = 1;
            break
        end


    end

    if failFlag ==0
        fprintf('Saving: %s\n',fullfile(L1Apath,strrep(fileList(i).name,'.TXT','L1A.mat')))
        save(fullfile(L1Apath,strrep(fileList(i).name,'.TXT','L1A.mat')), 'L1A')
        
        if makePlots
            plotFun(sensorList,L1A,strrep(fileList(i).name,'.TXT',''), flags)
        end
    end
    clear L1A badGPS badQ badRelAz badSZA badTilt T
end

%% Figures
function plotFun(sensorList,L1A,fileName, flags)
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

    exportgraphics(fh,sprintf('plt/%s_L1A_%s.png',fileName,(sensorList{s})))
    close
end
end





