% Process L1A DALEC to L1B
% 
% Apply calibrations and dark offsets
% Calibration file usage:
%
%    K1=d0*(V-DC)+d1
%    K2=e0*(V-DC)+e1
%    K3=f0*(V-DC)+f1
%    Ed=a0*((V-DC)/(Inttime+DeltaT_Ed)/K1)/(Tempco_Ed*(Temp-Tref)+1)
%    Lu=b0*((V-DC)/(Inttime+DeltaT_Lu)/K2)/(Tempco_Lu*(Temp-Tref)+1)
%    Lsky=c0*((V-DC)/(Inttime+DeltaT_Lsky)/K3)/(Tempco_Lsky*(Temp-Tref)+1)
%
% Interpolate to common timestamps and wavebands
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Inputs: L1A files from dalec_L1A.m, calibration file from IMO
%   Output: L1B files matlab structures
%           Combined plots of calibrated Ed, Lu, Lsky mean +/- std
%
% D. Aurin NASA/GSFC November 2024

% path(path,'./sub')            <-- uncomment if you're not me.
%% Setup
wipe
calPath = ...
    '~/Projects/HyperPACE/Instrument_Data_Cals/DALEC/DALEC_0012_2024_05_07_V3b.cal';% <-- Set this
dataPath = '~/Projects/HyperPACE/field_data/DALEC/NF2405_VIIRS';       % <-- Set this
L1Apath = fullfile(dataPath,'L1A/');
L1Bpath = fullfile(dataPath,'L1B/');
plotPath = fullfile(dataPath,'Plots/L1B/');

makePlots = 1;

fileList = dir(fullfile(L1Apath,'*.mat'));
sensorList = {'Ed','Lsky','Lu'};

if ~isfolder(L1Bpath)
    mkdir(L1Bpath);
end
% End setup
%% Read calibrations
% Read header:
fid = fopen(calPath,"r");
for l=1:50
    tline = fgetl(fid);
    if contains(tline,"; ")
        continue
    end
    if contains(tline," = ")
        disp(tline)
        keyValue = strsplit(tline,' = ');
        eval(sprintf('%s = %e;',keyValue{1},str2double(keyValue{2})))
    end
end
fclose(fid);

% Read Cals:
opts = detectImportOptions(calPath,...
    "FileType","text",...
    "Delimiter",",");
T = readtable(calPath,opts);

%% Read and process L1A 
for i=1:length(fileList)
    fpf = fullfile(fileList(i).folder,fileList(i).name);
    load(fpf) %L1A
    fprintf('Processing: %s\n',fpf)
    
    for s = 1:length(sensorList)
        fprintf('Applying sensor cals %s\n',sensorList{s})
        V = L1A.(sensorList{s}).light;
        DC = L1A.(sensorList{s}).dark;
        Inttime = L1A.(sensorList{s}).intTime;
        Temp = L1A.(sensorList{s}).temp;
        
        K1=d1*(V-DC)+d0;
        K2=e1*(V-DC)+e0;
        K3=f1*(V-DC)+f0;

        if strcmpi(sensorList{s},'ed')
            % (W/m^2/nm)
            Ed=T.a0' .* ((V-DC) ./ (Inttime+Delta_T_Ed) ./ K1)...
                ./ (T.Tempco_Ed' .* (Temp-Tref)+1);
        elseif strcmpi(sensorList{s},'lu')
            % (W/m^2/sr/nm)
            Lu=T.b0' .* ((V-DC) ./ (Inttime+Delta_T_Lu) ./ K2)...
                ./ (T.Tempco_Lu' .* (Temp-Tref)+1);
        elseif strcmpi(sensorList{s},'lsky')
            % (W/m^2/sr/nm)
            Lsky=T.c0' .* ((V-DC) ./ (Inttime+Delta_T_Lsky) ./ K3)...
                ./ (T.Tempco_Lsky' .* (Temp-Tref)+1);
        end        
        clear V DC Inttime Temp K1 K2 K3
    end
    
    % Interpolate in time
    % Lu is the slowest. Interpolate to Lu
    disp('Interpolating sensors to common timestamps (Lu)')
    dateTime = L1A.Lu.datetime;
    Ed = interp1(L1A.Ed.datetime, Ed,dateTime);    
    Lsky= interp1(L1A.Lsky.datetime, Lsky,dateTime);
    % Truncate NaNs (no extrapolation allowed)
    [rowsEd,~] = find(isnan(Ed));
    [rowsLsky,~] = find(isnan(Lsky));
    rowsBadTime = unique([rowsEd;rowsLsky]);
    if ~isempty(rowsBadTime)
        fprintf('Truncating %d timestamps\n',length(rowsBadTime))

        Lu(rowsBadTime,:) = [];
        Ed(rowsBadTime,:) = [];
        Lsky(rowsBadTime,:) = [];
        L1A.Lu.datetime(rowsBadTime) = [];
        L1A.Lu.lat(rowsBadTime) = [];
        L1A.Lu.lon(rowsBadTime) = [];
        L1A.Lu.sza(rowsBadTime) = [];
        L1A.Lu.relAz(rowsBadTime) = [];
        L1A.Lu.pitch(rowsBadTime) = [];
        L1A.Lu.roll(rowsBadTime) = [];
    end
    dateTime = L1A.Lu.datetime;

    % Interpolate to common wavebands
    % These appear equally spaced at ~3.3 nm
    disp('Interpolating spectrally to 3.3 nm')
    % wavelength = 340:3.3:964; % 190 elements
    wavelength = 340:3.3:932; % 190 elements
    Ed = interp1(T.Lambda_Ed, Ed', wavelength)';
    Lu = interp1(T.Lambda_Lu, Lu', wavelength)';
    Lsky = interp1(T.Lambda_Lsky, Lsky', wavelength)';
    % Check for out-of-range wavebands
    [~,colsLu] = find(isnan(Lu));
    [~,colsEd] = find(isnan(Ed));
    [~,colsLsky] = find(isnan(Lsky));
    colsBadWave = unique([colsLu;colsEd;colsLsky]);
    if ~isempty(colsBadWave)
        fprintf('%d columns found with NaNs. Check calibration and wavelength.\n',length(colsBadWave))
        break
    end

    % Build L1B
    L1B.datetime = dateTime;
    L1B.lat = L1A.Lu.lat;
    L1B.lon = L1A.Lu.lon;
    L1B.wavelength = wavelength;
    L1B.sza = L1A.Lu.sza;
    L1B.relAz = L1A.Lu.relAz;
    L1B.tilt = max(L1A.Lu.pitch, L1A.Lu.roll);
    L1B.Ed = Ed;
    L1B.Lu = Lu;
    L1B.Lsky = Lsky;
    clear L1A Ed Lu Lsky wavelength dateTime

    save(fullfile(L1Bpath,strrep(fileList(i).name,'L1A','L1B')), 'L1B')
    %% Figures
    if makePlots
        fileName = strrep(fileList(i).name,'L1A.mat','');
        fh = figure;
        for s =1:length(sensorList)
            subplot(3,1,s)
            plot(L1B.wavelength,mean(L1B.(sensorList{s})),'-k')
            hold on
            plot(L1B.wavelength, mean(L1B.(sensorList{s})) + std(L1B.(sensorList{s})),'--k')
            plot(L1B.wavelength,mean(L1B.(sensorList{s})) - std(L1B.(sensorList{s})),'--k')            
            grid on
            if strcmpi(sensorList{s},'ed')
                title(sprintf('%s: %s',sensorList{s},strrep(fileName,'_',' ')))
                ylabel('[W/m^2/nm]')
            else
                title(sensorList{s})
                ylabel('[W/m^2/sr/nm]')
            end
        end
        exportgraphics(fh,sprintf('%s%s_L1B.png',plotPath,fileName))
        close
    end    
    clear L1B
end