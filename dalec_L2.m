% Process L1B DALEC to L2
%
% Add ancillary data from field log file
% Break into hourly file groups for output
% Run the spectral filter (again, now that it's one-hour files)
% Extract 300s ensembles
% Drop brightest 90% of Lt(780)
% Take the slice mean of the ensemble for Lt,Li,Es
% Calculate the Zhang et al. 2017 glint correction
% Calculate Lw and Rrs for the ensemble with uncertainty
% Calculate the NIR residual correction
%   Only implemented for flat offset
%   TO DO: Improve by switching between flat offset and SimSpec based on AVW
% Export L2 ensembles as hourly files (mat) and SeaBASS files

wipe
%% Setup
ancPath = '~/Projects/HyperPACE/field_data/metadata/VIIRS2024/VIIRS2024_Ancillary.mat';
dataPath = '~/Projects/HyperPACE/field_data/DALEC/VIIRS2024';
L1Bpath = fullfile(dataPath,'L1B/');
L2path = fullfile(dataPath,'L2/');
sb.experiment = 'VIIRS_VALIDATION';
sb.cruise = 'NF2405_VIIRS';
sb.calibration_date = '20240507';
sb.L2path = L2path;

makePlots = 1;
makeSeaBASS = 1;

% Ancillary file matching
tLim = 15; % minutes (ship data interpolate from 10 to 5 minutes
dLim = 5; % km for stations. Field log lat/lon appear unreliable; increase threshold
kpdlat = 111.325; %km per deg. latitude

% Spectral filter
sigmaLight.Ed = 5;
sigmaLight.Lsky = 8;
sigmaLight.Lu = 3;
fRange = [350, 900];    % Spectral range for spectral filter

% Lowest X% of Lt (Hooker & Morel 2003; Hooker et al. 2002; Zibordi et al. 2002, IOCCG Protocols)
% depends on FOV and integration time of instrument. Hooker cites a rate of 2 Hz.
LtPrct = 10;

% Environment fallbacks when not present in ancillary data. All others are
% required in the L1B file
envDefaults.wind = 10;% wind[m/s]
envDefaults.aod = 0.15; % AOD(550)
envDefaults.cloud = 0; % cloud[%]
envDefaults.sst = 25; % sst[C]
envDefaults.sal = 34; % salt[PSU]

% Ensemble duration
ensSeconds = 300;

% Glint uncertainty
rhoUnc = 0.0017; %Estimated from FICE22

% NIR correction
applyNIR = 0;
NIRWave = [700 800];

% Negative Rrs range
testNegRrs = [400 700];

%% End setup
fileList = dir(fullfile(L1Bpath,'*.mat'));
sensorList = {'Ed','Lsky','Lu'};

if ~isfolder(L2path)
    mkdir(L2path);
end

%% Ancillary data
% No SST, AOD550, or salinity for this cruise. Only field notes.
load(ancPath) % ancillary ancHeaders

%% Read and process L1B

for i=1:length(fileList)
    fpf = fullfile(fileList(i).folder,fileList(i).name);
    load(fpf) %L1B
    fprintf('Processing: %s\n',fpf)

    %% Fold in ancillary data
    L1B.ancillary(length(L1B.datetime)) = struct();
    dateTime = num2cell(L1B.datetime);
    lat = num2cell(L1B.lat);
    lon = num2cell(L1B.lon);
    relAz = num2cell(L1B.relAz);
    sza = num2cell(L1B.sza);
    tilt = num2cell(L1B.tilt);
    [L1B.ancillary.datetime] = dateTime{:};
    [L1B.ancillary.lat] = lat{:};
    [L1B.ancillary.lon] = lon{:};
    [L1B.ancillary.relAz] = relAz{:};
    [L1B.ancillary.sza] = sza{:};
    [L1B.ancillary.tilt] = tilt{:};
    L1B = rmfield(L1B,{'datetime','lat','lon','relAz','sza','tilt'}); clear dateTime lat lon relAz sza tilt

    %% Break into hourly sections
    startHour = dateshift(L1B.ancillary(1).datetime,'start','hour');
    endHour = dateshift(L1B.ancillary(end).datetime,'end','hour');
    hourSteps = startHour:hours(1):endHour;

    for h=1:length(hourSteps)
        if h<length(hourSteps)
            match = find([L1B.ancillary.datetime] >= hourSteps(h) & ...
                [L1B.ancillary.datetime] < hourSteps(h+1));
        else
            match = find([L1B.ancillary.datetime] >= hourSteps(h));
        end

        if ~isempty(match)
            for n=match

                %% Match each L1B record to ancillary data if possible
                [nearTime, index] = find_nearest(L1B.ancillary(n).datetime,ancillary.datetime);
                % Within tLim min and stationary
                tDiff = abs(L1B.ancillary(n).datetime - nearTime);
                if tDiff <= minutes(tLim)
                    kpdlon = kpdlat*cos(pi*L1B.ancillary(n).lat/180);
                    dlat = kpdlat*(L1B.ancillary(n).lat - ancillary.lat(index));
                    dlon = kpdlon*(L1B.ancillary(n).lon - ancillary.lon(index));
                    dist = sqrt(dlat.^2 + dlon.^2); % distance to station [km]
                    if dist < dLim
                        flagMatch = 1;
                        L1B.ancillary(n).station = ancillary.station(index);
                        L1B.ancillary(n).wind = ancillary.wind(index);
                        L1B.ancillary(n).cloud = ancillary.cloud(index);
                        L1B.ancillary(n).waves = ancillary.waves(index);
                        % fprintf('Match Station: %d\n', ancillary.station(index))
                        L1B.ancillary(n).sst = envDefaults.sst; % Not present in ancillary
                        L1B.ancillary(n).sal = envDefaults.sal; % Not present in ancillary
                        L1B.ancillary(n).aod = envDefaults.aod; % Not present in ancillary
                    else
                        flagMatch = 0;
                    end
                else
                    flagMatch = 0;
                end
                if ~flagMatch
                    L1B.ancillary(n).station = nan;
                    L1B.ancillary(n).wind = envDefaults.wind;
                    L1B.ancillary(n).cloud = envDefaults.cloud;
                    L1B.ancillary(n).waves = nan;
                    L1B.ancillary(n).sst = envDefaults.sst;
                    L1B.ancillary(n).sal = envDefaults.sal;
                    L1B.ancillary(n).aod = envDefaults.aod;
                end
            end
            %% Build hourly structure based on match
            Hourly.ancillary(length(match)) = struct();
            [Hourly.ancillary.datetime] = L1B.ancillary(match).datetime;
            [Hourly.ancillary.lat] = L1B.ancillary(match).lat;
            [Hourly.ancillary.lon] = L1B.ancillary(match).lon;
            [Hourly.ancillary.relAz] = L1B.ancillary(match).relAz;
            [Hourly.ancillary.sza] = L1B.ancillary(match).sza;
            [Hourly.ancillary.tilt] = L1B.ancillary(match).tilt;
            [Hourly.ancillary.station] = L1B.ancillary(match).station;
            [Hourly.ancillary.wind] = L1B.ancillary(match).wind;
            [Hourly.ancillary.cloud] = L1B.ancillary(match).cloud;
            [Hourly.ancillary.waves] = L1B.ancillary(match).waves;
            [Hourly.ancillary.sst] = L1B.ancillary(match).sst;
            [Hourly.ancillary.sal] = L1B.ancillary(match).sal;
            [Hourly.ancillary.aod] = L1B.ancillary(match).aod;
            [Hourly.Ed] = L1B.Ed(match,:);
            [Hourly.Lu] = L1B.Lu(match,:);
            [Hourly.Lsky] = L1B.Lsky(match,:);
            clear match

            %% Run spectral filter on hourly file
            for s = 1:length(sensorList)
                [~,wv1] = find_nearest(fRange(1),L1B.wavelength);
                [~,wv2] = find_nearest(fRange(2),L1B.wavelength);
                subLight = Hourly.(sensorList{s})(:,wv1:wv2);
                peakLight = max(Hourly.(sensorList{s})(:,wv1:wv2),[],2);
                normLight = subLight ./ peakLight;

                ave = mean(normLight);
                stD = std(normLight);
                sigma = sigmaLight.(sensorList{s});
                flag = normLight < (ave - sigma*stD) |...
                    normLight > (ave + sigma*stD);
                [badSpec,~] = find(flag);
                badSpec = unique(badSpec);
                fprintf('%s records of %d removed for spectral outliers: %d\n', sensorList{s}, length(flag), length(badSpec))
                Hourly.Ed(badSpec,:) = [];
                Hourly.Lu(badSpec,:) = [];
                Hourly.Lsky(badSpec,:) = [];
                Hourly.ancillary(badSpec) = [];
                clear normLight ave stD flag
            end

            %% Extract ensembles
            startEns = dateshift(Hourly.ancillary(1).datetime,'start','minute');
            endEns = dateshift(Hourly.ancillary(end).datetime,'end','minute');
            ensSteps = startEns:seconds(ensSeconds):endEns;

            ens = 0;
            for e=1:length(ensSteps)
                if e<length(ensSteps)
                    match = find([Hourly.ancillary.datetime] >= ensSteps(e) & ...
                        [Hourly.ancillary.datetime] < ensSteps(e+1));
                else
                    match = find([Hourly.ancillary.datetime] >= ensSteps(e));
                end

                if ~isempty(match)
                    ens = ens+1;

                    %% Glitter correction
                    Lt = Hourly.Lu(match,:);
                    [wv780,whr780] = find_nearest(780,L1B.wavelength);
                    Lt780 = Lt(:,whr780);
                    % There have to be at least 20 to get 2 left at 10% Lt
                    if length(Lt780>=20)
                        [sortLt,indxLt] = sort(Lt780);
                        x = round(length(Lt780) * LtPrct / 100);
                        xIndx = indxLt(1:x);
                        % Redefine the match to only include lowest 10% Lt
                        match = match(xIndx);
                        N = length(match);
                        fprintf('Number of ensemble records after glitter reduction: %d\n',N)

                        L2.Radiance.Lt(ens,:) = mean(Hourly.Lu(match,:)); % Change Lu to Lt
                        L2.Irradiance.Es(ens,:) = mean(Hourly.Ed(match,:)); % Change Ed to Es
                        L2.Radiance.Li(ens,:) = mean(Hourly.Lsky(match,:)); % Change Lsky to Li
                        L2.Radiance.Lt_sd(ens,:) = std(Hourly.Lu(match,:));
                        L2.Irradiance.Es_sd(ens,:) = std(Hourly.Ed(match,:));
                        L2.Radiance.Li_sd(ens,:) = std(Hourly.Lsky(match,:));

                        L2.wavelength = L1B.wavelength;

                        L2.ancillary(ens).datetime = mean([Hourly.ancillary(match).datetime]);
                        L2.ancillary(ens).lat = mean([Hourly.ancillary(match).lat]);
                        L2.ancillary(ens).lon = mean([Hourly.ancillary(match).lon]);
                        L2.ancillary(ens).relAz = mean([Hourly.ancillary(match).relAz]);
                        L2.ancillary(ens).sza = mean([Hourly.ancillary(match).sza]);
                        L2.ancillary(ens).tilt = mean([Hourly.ancillary(match).tilt]);

                        L2.ancillary(ens).station = nanmean([Hourly.ancillary(match).station]);
                        L2.ancillary(ens).wind = nanmean([Hourly.ancillary(match).wind]);
                        L2.ancillary(ens).cloud = nanmean([Hourly.ancillary(match).cloud]);
                        L2.ancillary(ens).waves = nanmean([Hourly.ancillary(match).waves]);
                        L2.ancillary(ens).sst = nanmean([Hourly.ancillary(match).sst]);
                        L2.ancillary(ens).sal = nanmean([Hourly.ancillary(match).sal]);
                        L2.ancillary(ens).aod = nanmean([Hourly.ancillary(match).aod]);
                        L2.ancillary(ens).binCount = N;

                        %% Glint Correction
                        % Zhang, X., S. He, A. Shabani, P.-W. Zhai, and K. Du. 2017. Spectral sea
                        % surface reflectance of skylight. Opt. Express 25: A1-A13,
                        % doi:10.1364/OE.25.0000A1.
                        %   NOTE: dat/db.mat is a link to a local database.
                        %   Others will need to copy or link from the database
                        %   in HyperCP to this location.
                        disp('Calculating Z17 rho')

                        % Truncate to model limits for wind, AOD, SZA
                        wind = L2.ancillary(ens).wind;
                        if wind > 15
                            wind = 15;
                        end
                        aod = L2.ancillary(ens).aod;
                        if aod > 0.2
                            aod = 0.2;
                        end
                        sza = L2.ancillary(ens).sza;
                        if sza > 60
                            sza = 60;
                        end

                        % === environmental conditions during experiment ===
                        env.wind = wind; % wind speeds in m/s 0:5:15
                        env.od = aod; % aersosol optical depth at 550 nm 0:0.05:0.2
                        env.C = L2.ancillary(ens).cloud; % cloud cover. 0 (not used)
                        env.zen_sun = sza; % sun zenith angle 0:10:60
                        env.wtem = L2.ancillary(ens).sst; % water temperature (Deg C)
                        env.sal = L2.ancillary(ens).sal; % salinity PSU

                        % === The sensor ===
                        % the zenith and azimuth angles of light that the sensor will see
                        % 0 azimuth angle is where the sun located
                        % positive z is upward
                        % Relative azimuth convention: sensor to 180-sun azimuth
                        relAz = 180 - L2.ancillary(ens).relAz;
                        sensor.ang = [40, relAz]; % zenith and azimuth angle in degree
                        sensor.wv = L1B.wavelength; % wavelength in nm 250:5:1000
                        sensor.ang2 = sensor.ang + [0, 180]; % location where skylight is measured

                        rho = get_sky_sun_rho(env, sensor);
                        % Screen for NaNs (out of model range)
                        badRho = isnan(rho.rho);

                        rho_unc = repmat(rhoUnc,1,length(L2.wavelength)); %Estimated at 0.0017 from FICE22

                        Lw =  L2.Radiance.Lt(ens,:) - (rho.rho .* L2.Radiance.Li(ens,:));
                        % Lw uncertainty
                        Lw_unc = sqrt( (L2.Radiance.Lt_sd(ens,:)).^2 +...
                            (L2.Radiance.Li_sd(ens,:)).^2);% + ...
                        % (rho_unc).^2);

                        Rrs = Lw ./ L2.Irradiance.Es(ens,:);
                        % Rrs uncertainty
                        Rrs_unc = abs(Rrs) .* ...
                            sqrt( (Lw_unc./Lw).^2 +...
                            (L2.Irradiance.Es_sd(ens,:)./L2.Irradiance.Es(ens,:)).^2 );

                        if sum(badRho) < 0.5*length(badRho)
                            Lw(badRho) = [];
                            Lw_unc(badRho) = [];
                            Rrs(badRho) = [];
                            Rrs_unc(badRho) = [];


                            L2.Radiance.Lw(ens,:) = Lw;
                            L2.Radiance.Lw_unc(ens,:) = Lw_unc;
                            L2.Reflectance.Rrs(ens,:) = Rrs;
                            L2.Reflectance.Rrs_unc(ens,:) = Rrs_unc;

                            L2.wavelength_rho = L2.wavelength(~badRho);
                            % plot(L2.wavelength_rho,L2.Reflectance.Rrs_unc)

                            [L2.Reflectance.QWIP(ens),QCI,L2.Reflectance.AVW(ens)] = AVW_QWIP_2D_fun(...
                                L2.Reflectance.Rrs(ens,:),L2.wavelength_rho,'none','none');
                            %% NIR Residual Correction
                            if applyNIR
                                % Standard
                                % The approach in HyperCP was adapted from taking a NIR
                                % average offset to a minimum NIR value 700-800 nm
                                whrWave = L2.wavelength > NIRWave(1) & L2.wavelength <= NIRWave(2);
                                NIR = min(L2.Reflectance.Rrs(ens,whrWave));

                                L2.Reflectance.Rrs(ens,:) = L2.Reflectance.Rrs(ens,:) - NIR;

                                % SimSpec (PLACEHOLDER)
                                rho = L2.Reflectance.Rrs_unc(ens,:) * pi;
                            end

                            %% QC for negative Rrs (400:700)
                            whrWave = L2.wavelength > testNegRrs(1) & L2.wavelength <= testNegRrs(2);
                            negVIS = L2.Reflectance.Rrs(ens,whrWave) < 0;
                            if any(negVIS)
                                fprintf('Negative spectra removed ensemble %d\n',ens)
                                % L2.ancillary(ens) = [];
                                % L2.Reflectance.Rrs(ens,:) = [];
                                % L2.Reflectance.Rrs_unc(ens,:) = [];
                                ens = ens-1;
                            end
                        else
                            fprintf('Zhang glint model failure ensemble %d\n',ens)
                            ens = ens -1;
                            continue
                        end
                    else
                        fprintf('Too few spectra remaining ensemble %d\n',ens)
                        ens = ens -1;
                        continue
                    end

                end
            end

            %% Save hourly file of ensembles
            firstEns = datetime(L2.ancillary(1).datetime,'format','yyyyMMdd''T''HHmmss');
            outFile.L1B = fileList(i).name;
            outFile.L2 = sprintf('%s%s_L2.mat',fileList(i).name(1:9),string(firstEns));
            fprintf('Saving hourly file: %s\n',outFile.L2)
            save(fullfile(L2path,outFile.L2),"L2")

            %% Figures
            if makePlots
                plotL2(L2, dataPath, outFile.L2)
                close all
            end

            %% Figures
            if makeSeaBASS
                writeSeaBASS(L2, outFile, sb)
            end
            clear Hourly L2
        end

    end
end
