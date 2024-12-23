% Process L1B DALEC to L2
%
% Add ancillary data from field log file       
%   You need to build an ancillary structure in a file with:  <--- Note!!
%       Required fields: datetime, lat, lon,
%       Strongly encouraged: station, wind, sst, sal, AOD, cloud, waves
%   This will be matched to your radiometry as described below
%
% Break into hourly file groups for output
% Run the spectral outlier filter
% Extract 300s ensembles
% Drop brightest 90% of Lt(780) for glitter
% Take the slice mean of the ensemble for Lt,Li,Es
%
% Calculate the Zhang et al. 2017 glint correction
%       This will require the database which can be downloaded here: <--- Note!!
%           https://oceancolor.gsfc.nasa.gov/fileshare/dirk_aurin/
%       It's ~2 GB and needs to be renamed and placed in ./dat/db.mat
%
% Calculate Lw and Rrs for the ensemble with uncertainty
% 
% Calculate the NIR residual correction
%   Currently only implemented for flat offset
%   TO DO: Improve by switching between flat offset and SimSpec based on AVW
%
% Screen for negative Rrs
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Inputs: L1B files from dalec_L1B.m
%   Output: L2 files Hourly .mat files, L2 Hourly .sb files
%           Combined plots of calibrated Es, Lt, Lw, Rrs ensembles mean +/-
%           std with some metadata on QC
%
% D. Aurin, NASA/GSFC November 2024

% path(path,'./sub')                    <-- uncomment if you're not me.
%% Setup
wipe
% Ancillary data structure:
ancPath = 'dat/VIIRS2024_Ancillary.mat'; % <-- Set this
% ancPath = ...
%     '~/Projects/HyperPACE/field_data/metadata/VIIRS2024/VIIRS2024_Ancillary.mat'; % <-- Set this
dataPath = '~/Projects/HyperPACE/field_data/DALEC/NF2405_VIIRS';                    % <-- Set this
L1Bpath = fullfile(dataPath,'L1B/');
L2path = fullfile(dataPath,'L2/');

sb.experiment = 'VIIRS_VALIDATION';     % <-- Set this
sb.cruise = 'NF2405_VIIRS';             % <-- Set this
sb.calibration_date = '20240507';       % <-- Set this
sb.L2path = L2path;

makePlots = 1;
makeSeaBASS = 1;

% Ancillary file matching criteria
tLim = 15; % [minutes]
dLim = 5; % distance [km] for stations. Usually < 1 km. Field log lat/lon appear unreliable here; increase threshold
kpdlat = 111.325; %km per deg. latitude

% Spectral filter Sigmas
sigmaLight.Ed = 5;
sigmaLight.Lsky = 8;
sigmaLight.Lu = 3;
fRange = [350, 900];    % Spectral range for outlier filter

% Lowest X% of Lt (Hooker & Morel 2003; Hooker et al. 2002; Zibordi et al. 2002, IOCCG Protocols)
% depends on FOV and integration time of instrument. Hooker cites a rate of
% 2 Hz for multispectral instruments. Hyper is slower (generally < 0.5 Hz).
LtPrct = 10;

% Environment fallbacks when not present in ancillary data. All others are
% required in the L1B file
envDefaults.wind = 5;       % wind[m/s]            % <-- Set this
envDefaults.aod = 0.18;     % AOD(550)             % <-- Set this
envDefaults.cloud = nan;    % cloud[%]             % <-- Set this
envDefaults.sst = 27;       % sst[C]               % <-- Set this
envDefaults.sal = 36;       % salt[PSU]            % <-- Set this

% Ensemble duration
ensSeconds = 300;   % HyperCP default, but this depends on your sampling, speed, stations, etc.

% Glint uncertainty
% rhoUnc = 0.0017; %Estimated from FICE22
rhoUnc = 0.003; % Estimated from Ruddick et al. 2006

% NIR correction
applyNIR = 1 ;                                   % <-- Set this
NIRWave = [700 800];

% Negative Rrs spectral range to test
testNegRrs = [395 715]; %UVA to NIR

%% End setup
%% Process L2
if ~isfolder(L2path)
    mkdir(L2path);
end
fileList = dir(fullfile(L1Bpath,'*.mat'));
sensorList = {'Ed','Lsky','Lu'};

%% Ancillary data
% Currently no SST, AOD550, or salinity for this cruise. Only field notes.
load(ancPath) % ancillary ancHeaders
% Pack with nans where necessary
ancOptFields = {'station','wind', 'sst', 'sal', 'aod', 'cloud', 'waves'};
for i=1:length(ancOptFields)
    if ~isfield(ancillary,ancOptFields{i})
        ancillary.(ancOptFields{i}) = nan*ancillary.lat;
    end
end
% Now backfill with defaults (for matches)
ancillary.wind(isnan(ancillary.wind)) = envDefaults.wind;
ancillary.sst(isnan(ancillary.sst)) = envDefaults.sst;
ancillary.sal(isnan(ancillary.sal)) = envDefaults.sal;
ancillary.aod(isnan(ancillary.aod)) = envDefaults.aod;
ancillary.cloud(isnan(ancillary.cloud)) = envDefaults.cloud;

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
            matchHour = find([L1B.ancillary.datetime] >= hourSteps(h) & ...
                [L1B.ancillary.datetime] < hourSteps(h+1));
        else
            matchHour = find([L1B.ancillary.datetime] >= hourSteps(h));
        end

        if ~isempty(matchHour)
            for n=matchHour
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
                        L1B.ancillary(n).sst = ancillary.sst(index);
                        L1B.ancillary(n).sal = ancillary.sal(index);
                        L1B.ancillary(n).aod = ancillary.aod(index);
                    else
                        flagMatch = 0;
                    end
                else
                    flagMatch = 0;
                end
                if ~flagMatch
                    L1B.ancillary(n).station = nan; % Not required for processing                    
                    L1B.ancillary(n).cloud = nan;% Not required for processing
                    L1B.ancillary(n).waves = nan; % Not required for processing
                    L1B.ancillary(n).wind = envDefaults.wind;
                    L1B.ancillary(n).sst = envDefaults.sst;
                    L1B.ancillary(n).sal = envDefaults.sal;
                    L1B.ancillary(n).aod = envDefaults.aod;
                end
            end            

            %% Build hourly structure based on matchHour
            Hourly.ancillary(length(matchHour)) = struct();
            [Hourly.ancillary.datetime] = L1B.ancillary(matchHour).datetime;
            [Hourly.ancillary.lat] = L1B.ancillary(matchHour).lat;
            [Hourly.ancillary.lon] = L1B.ancillary(matchHour).lon;
            [Hourly.ancillary.relAz] = L1B.ancillary(matchHour).relAz;
            [Hourly.ancillary.sza] = L1B.ancillary(matchHour).sza;
            [Hourly.ancillary.tilt] = L1B.ancillary(matchHour).tilt;
            [Hourly.ancillary.station] = L1B.ancillary(matchHour).station;
            [Hourly.ancillary.wind] = L1B.ancillary(matchHour).wind;
            [Hourly.ancillary.cloud] = L1B.ancillary(matchHour).cloud;
            [Hourly.ancillary.waves] = L1B.ancillary(matchHour).waves;
            [Hourly.ancillary.sst] = L1B.ancillary(matchHour).sst;
            [Hourly.ancillary.sal] = L1B.ancillary(matchHour).sal;
            [Hourly.ancillary.aod] = L1B.ancillary(matchHour).aod;
            [Hourly.Ed] = L1B.Ed(matchHour,:);
            [Hourly.Lu] = L1B.Lu(matchHour,:);
            [Hourly.Lsky] = L1B.Lsky(matchHour,:);
            clear matchHour

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
                Hourly.specFilter.(sensorList{s}).startN = length(flag); % May be reduced with each sensor
                Hourly.specFilter.(sensorList{s}).Nremoved = length(badSpec);
                fprintf('Hourly %s records removed for spectral outliers: %d / %d\n', sensorList{s}, length(badSpec), length(flag))
                Hourly.Ed(badSpec,:) = [];
                Hourly.Lu(badSpec,:) = [];
                Hourly.Lsky(badSpec,:) = [];
                Hourly.ancillary(badSpec) = [];
                clear normLight ave stD flag
            end

            %% Extract ensembles
            emptyFlag = 0;  % In case all ensembles are lost
            nRrs = 0;       % Track ensembles removed for negative Rrs
            L2.QC.negRrsRemoved = nRrs; % This will increment as needed
            startEns = dateshift(Hourly.ancillary(1).datetime,'start','minute');
            endEns = dateshift(Hourly.ancillary(end).datetime,'end','minute');
            ensSteps = startEns:seconds(ensSeconds):endEns;

            ens = 0;
            for e=1:length(ensSteps)
                if e<length(ensSteps)
                    matchEns = find([Hourly.ancillary.datetime] >= ensSteps(e) & ...
                        [Hourly.ancillary.datetime] < ensSteps(e+1));
                else
                    matchEns = find([Hourly.ancillary.datetime] >= ensSteps(e));
                end

                if ~isempty(matchEns)   
                    fprintf('Ensemble spectra initially: %d\n',length(matchEns))
                    %% Glitter correction to the ensemble
                    Lt = Hourly.Lu(matchEns,:);
                    [wv780,whr780] = find_nearest(780,L1B.wavelength);
                    Lt780 = Lt(:,whr780);

                    % There have to be at least 20 to get 2 left at 10% Lt
                    if length(Lt780) >= 20
                        
                        ens = ens+1;
                        
                        [sortLt,indxLt] = sort(Lt780);
                        x = round(length(Lt780) * LtPrct / 100);
                        xIndx = indxLt(1:x);
                        % Redefine the matchEns to only include lowest 10% Lt
                        matchEns = matchEns(xIndx);
                        N = length(matchEns);
                        fprintf('Number of ensemble records after glitter reduction: %d\n',N)
                        L2.QC.glitterStats(ens).startN = length(Lt780);
                        L2.QC.glitterStats(ens).finalN = N;

                        L2.Radiance.Lt(ens,:) = mean(Hourly.Lu(matchEns,:)); % Change Lu to Lt
                        L2.Irradiance.Es(ens,:) = mean(Hourly.Ed(matchEns,:)); % Change Ed to Es
                        L2.Radiance.Li(ens,:) = mean(Hourly.Lsky(matchEns,:)); % Change Lsky to Li
                        L2.Radiance.Lt_sd(ens,:) = std(Hourly.Lu(matchEns,:));
                        L2.Irradiance.Es_sd(ens,:) = std(Hourly.Ed(matchEns,:));
                        L2.Radiance.Li_sd(ens,:) = std(Hourly.Lsky(matchEns,:));

                        if ens==1
                            L2.QC.specFilter.Lt_Hourly_Nstart = Hourly.specFilter.Lu.startN;
                            L2.QC.specFilter.Lt_Hourly_Nremoved = Hourly.specFilter.Lu.Nremoved;
                            L2.QC.specFilter.Li_Hourly_Nstart = Hourly.specFilter.Lsky.startN;
                            L2.QC.specFilter.Li_Hourly_Nremoved = Hourly.specFilter.Lsky.Nremoved;
                            L2.QC.specFilter.Es_Hourly_Nstart = Hourly.specFilter.Ed.startN;
                            L2.QC.specFilter.Es_Hourly_Nremoved = Hourly.specFilter.Ed.Nremoved;
                        end

                        L2.wavelength = L1B.wavelength;

                        L2.ancillary(ens).rawFileName = strrep(fileList(i).name,'L1B.mat','.TXT');
                        L2.ancillary(ens).datetime = mean([Hourly.ancillary(matchEns).datetime]);
                        L2.ancillary(ens).lat = mean([Hourly.ancillary(matchEns).lat]);
                        L2.ancillary(ens).lon = mean([Hourly.ancillary(matchEns).lon]);
                        L2.ancillary(ens).relAz = mean([Hourly.ancillary(matchEns).relAz]);
                        L2.ancillary(ens).sza = mean([Hourly.ancillary(matchEns).sza]);
                        L2.ancillary(ens).tilt = mean([Hourly.ancillary(matchEns).tilt]);

                        L2.ancillary(ens).station = nanmean([Hourly.ancillary(matchEns).station]);
                        L2.ancillary(ens).wind = nanmean([Hourly.ancillary(matchEns).wind]);
                        L2.ancillary(ens).cloud = nanmean([Hourly.ancillary(matchEns).cloud]);
                        L2.ancillary(ens).waves = nanmean([Hourly.ancillary(matchEns).waves]);
                        L2.ancillary(ens).sst = nanmean([Hourly.ancillary(matchEns).sst]);
                        L2.ancillary(ens).sal = nanmean([Hourly.ancillary(matchEns).sal]);
                        L2.ancillary(ens).aod = nanmean([Hourly.ancillary(matchEns).aod]);
                        L2.ancillary(ens).binCount = N;

                        %% Glint Correction
                        % Zhang, X., S. He, A. Shabani, P.-W. Zhai, and K. Du. 2017. Spectral sea
                        % surface reflectance of skylight. Opt. Express 25: A1-A13,
                        % doi:10.1364/OE.25.0000A1.
                        %   NOTE: dat/db.mat is a link to a local database.
                        %   Others will need to copy or link from the database
                        %   in HyperCP to this location.
                        fprintf('Calculating Z17 rho for ensemble %d\n',ens)

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
                        env.C = L2.ancillary(ens).cloud; % cloud cover. (not used)
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
                        nanRho = isnan(rho.rho);

                        rho_unc = repmat(rhoUnc,1,length(L2.wavelength)); %Estimated at 0.0017 from FICE22

                        Lw =  L2.Radiance.Lt(ens,:) - (rho.rho .* L2.Radiance.Li(ens,:));
                        % Lw uncertainty (assuming random, uncorrelated error)
                        Lw_unc =  ...
                            sqrt(...
                            ( (L2.Radiance.Lt_sd(ens,:)).^2 +...
                            (rho_unc.*L2.Radiance.Li(ens,:)).^2 + ...
                            (rho.rho.*L2.Radiance.Li_sd(ens,:)).^2 )...
                            );

                        Rrs = Lw ./ L2.Irradiance.Es(ens,:);
                        % Rrs uncertainty
                        % Rrs_unc = abs(Rrs) .* ...
                        %     sqrt(...
                        %     (Lw_unc./Lw).^2 + ...
                        %     (L2.Irradiance.Es_sd(ens,:)./L2.Irradiance.Es(ens,:)).^2 );
                        %% The sensitivity coefficients are missing here and this needs to be revisited.
                        Rrs_unc = abs(Rrs) .* ...
                            sqrt(...
                            (L2.Radiance.Li_sd(ens,:)./L2.Radiance.Li(ens,:)).^2 + ...
                            (rho_unc./rho.rho).^2 + ...
                            (L2.Radiance.Lt_sd(ens,:)./L2.Radiance.Lt(ens,:)).^2 + ...
                            (L2.Irradiance.Es_sd(ens,:)./L2.Irradiance.Es(ens,:)).^2 );

                        % Only retain if more than 80% of rho survived the
                        % Zhang model
                        if sum(nanRho) < 0.2*length(nanRho)

                            % Eliminate NaN spectral range for Lw and Rrs
                            Lw(nanRho) = [];
                            Lw_unc(nanRho) = [];
                            rho.rho(nanRho) = [];
                            Rrs(nanRho) = [];
                            Rrs_unc(nanRho) = [];                            

                            L2.Radiance.Lw(ens,:) = Lw;
                            L2.Radiance.Lw_unc(ens,:) = Lw_unc;
                            L2.Radiance.rho(ens,:) = rho.rho;
                            L2.Radiance.rho_unc(ens,:) = rhoUnc;
                            L2.Reflectance.Rrs(ens,:) = Rrs;
                            L2.Reflectance.Rrs_unc(ens,:) = Rrs_unc;

                            L2.wavelength_rho = L2.wavelength(~nanRho);
                            % plot(L2.wavelength_rho,L2.Reflectance.Rrs_unc)

                            [L2.Reflectance.QWIP(ens),QCI,L2.Reflectance.AVW(ens)] = AVW_QWIP_2D_fun(...
                                L2.Reflectance.Rrs(ens,:),L2.wavelength_rho,'none','none');

                            %% NIR Residual Correction
                            if applyNIR
                                % Standard flat NIR min value.
                                % The approach in HyperCP was adapted from taking a NIR
                                % average offset to a minimum NIR value
                                % 700-800 nm. Debatable...
                                whrWave = L2.wavelength > NIRWave(1) & L2.wavelength <= NIRWave(2);
                                L2.Reflectance.NIR(ens) = min(L2.Reflectance.Rrs(ens,whrWave));

                                L2.Reflectance.Rrs(ens,:) = L2.Reflectance.Rrs(ens,:) - L2.Reflectance.NIR(ens);

                                % SimSpec (PLACEHOLDER for optically turbid waters)
                                % rhos = L2.Reflectance.Rrs_unc(ens,:) * pi;
                            end                            

                            %% QC for negative Rrs in the VIS (400:700)
                            whrWave = L2.wavelength > testNegRrs(1) & L2.wavelength <= testNegRrs(2);
                            negVIS = L2.Reflectance.Rrs(ens,whrWave) < 0;
                            if any(negVIS)
                                fprintf('Negative spectra removed ensemble %d\n',ens)
                                nRrs = nRrs +1;
                                L2.QC.negRrsRemoved = nRrs;                                

                                % Need to terminate if this is final ensemble...
                                if e==length(ensSteps)
                                    [L2, emptyFlag] = deleteLastEns(L2, ens);
                                    break
                                end
                                ens = ens-1;
                            else
                                % Set remaining negative reflectance to NaN (UV and NIR only!)
                                L2.Reflectance.Rrs(ens,(L2.Reflectance.Rrs(ens,:) < 0)) = nan;
                            end

                            
                        else
                            fprintf('Zhang glint model failure ensemble %d\n',ens)
                            % Need to terminate if this is final ensemble...
                            if e==length(ensSteps)                                
                                [L2, emptyFlag] = deleteLastEns(L2, ens);
                                break
                            end
                            ens = ens -1;
                        end
                    else
                        fprintf('Too few spectra remaining: %d\n',length(Lt780))
                        % Need to terminate if this is final ensemble...
                        if e==length(ensSteps)
                            [L2, emptyFlag] = deleteLastEns(L2, ens);
                            break
                        end
                        ens = ens -1;
                    end

                else
                    % ens has not been incremented, nor L2 populated for ens at this point.
                    fprintf('No matching spectra ensemble: %d\n',ens+1)
                    % Need to terminate if this is final ensemble...
                    if e==length(ensSteps)
                        % [L2, emptyFlag] = deleteLastEns(L2, ens);
                        break
                    end
                    % ens = ens -1;
                end
            end

            %% Save hourly file of ensembles            
            if ~emptyFlag
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

                %% SeaBASS File Writer
                if makeSeaBASS
                    writeSeaBASS(L2, outFile, sb, {'Rrs','Es'})
                end

            else
                disp('FAIL: No ensembles remaining. No file saved.')
            end
            clear Hourly L2
        end

    end
end

%%
function [L2, emptyFlag] = deleteLastEns(L2, ens)

emptyFlag = 0;
if ens ==0
    disp('Not expecting this anymore ################################')
    % ens =1;
elseif ens == 1
    % This can happen if all ensembles failed QC and the last ensemble has no matching spectra
    emptyFlag = 1;
end
% Remove this final ensemble
L2.ancillary(ens) = [];
L2.QC.glitterStats(ens) = [];
% L2.QC.specFilter(ens) = []; This is hourly, not ensemble
groups = {'Radiance','Irradiance','Reflectance'};
for g=1:length(groups)
    fieldNames = fieldnames(L2.(groups{g}));
    for f=1:length(fieldNames)
        if min(size(L2.(groups{g}).(fieldNames{f}))) ==1
            % fprintf('%s %s\n',groups{g},fieldNames{f})
            L2.(groups{g}).(fieldNames{f})(ens) = [];
        else
            % fprintf('%s %s\n',groups{g},fieldNames{f})
            L2.(groups{g}).(fieldNames{f})(ens,:) = [];
        end
    end
end

end