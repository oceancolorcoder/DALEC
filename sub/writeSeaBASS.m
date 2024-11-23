function writeSeaBASS(L2, outFile, sb, dataTypes)

if ~isfolder(sprintf('%s/SeaBASS',sb.L2path))
    mkdir(sprintf('%s/SeaBASS',sb.L2path))
end

% Make a template seabass header for your data and place in /dat
inHeader = sprintf('dat/%s_L2_header.sb', sb.cruise);
[~, header] = readsb(inHeader);

header = rmfield(header,'wavelength_lists');
header = rmfield(header,'fields_list');
hFields = fields(header);


outFile.SB = strrep(outFile.L2,'DALEC_12_','');
outFile.SB = sprintf('%s_%s_DALEC_%s', sb.experiment, sb.cruise, strrep(outFile.SB,'.mat','.sb'));
header.data_file_name = outFile.SB;
outFile.SB = sprintf('%sSeaBASS/%s', sb.L2path, outFile.SB);
outFile.Raw = strrep(outFile.L1B,'L1B.mat','.raw');

nullChar = header.missing;
startTime = min([L2.ancillary.datetime]);
stopTime = max([L2.ancillary.datetime]);
southLat = min([L2.ancillary.lat]);
northLat = max([L2.ancillary.lat]);
westLon = min([L2.ancillary.lon]);
eastLon = max([L2.ancillary.lon]);

header.start_date = char(datetime(startTime,'Format','yyyyMMdd'));
header.start_date = char(datetime(startTime,'Format','yyyyMMdd'));
header.start_time = char(datetime(startTime,'Format','hh:mm:ss'));
header.end_date = char(datetime(stopTime,'Format','yyyyMMdd'));
header.end_time = char(datetime(stopTime,'Format','hh:mm:ss'));
header.north_latitude = northLat;
header.south_latitude = southLat;
header.west_longitude = westLon;
header.east_longitude = eastLon;

header.original_file_name = outFile.Raw;

% Fix nans
for i=1:length(hFields)
    if ~ischar(header.(hFields{i}))
        if isnan(header.(hFields{i}))
            header.(hFields{i}) = nullChar;
        end
    end
end

ancCols = 17; % Includes breakout of date(y M d) and time (h m s)

for dT = 1:length(dataTypes)
    if strcmpi(dataTypes{dT},'Rrs')
        SBfields = 'station,date,time,lat,lon,RelAz,SZA,cloud,wind,waveht,sal,sst,bincount';
        SBunits = 'none,yyyyMMdd,hh:mm:ss,degrees,degrees,degrees,degrees,%,m/s,m,psu,degreesC,none';
        %% Rrs
        wavelength = L2.wavelength_rho;
        for i=1:length(wavelength)
            SBfields = [SBfields  sprintf(',Rrs%.1f',wavelength(i))];
            SBunits = [SBunits ',1/sr'];
        end
        for i=1:length(wavelength)
            SBfields = [SBfields  sprintf(',Rrs%.1f_unc',wavelength(i))];
            SBunits = [SBunits ',1/sr'];
        end

        rows = length([L2.ancillary.datetime]);
        cols = ancCols + 2*length(wavelength);
        dataBlock = nan(rows,cols);
        for i=1:rows
            dataBlock(i,:) = [
                L2.ancillary(i).station year(L2.ancillary(i).datetime) month(L2.ancillary(i).datetime) day(L2.ancillary(i).datetime) ...
                hour(L2.ancillary(i).datetime) minute(L2.ancillary(i).datetime) round(second(L2.ancillary(i).datetime)) ...
                L2.ancillary(i).lat L2.ancillary(i).lon L2.ancillary(i).relAz  L2.ancillary(i).sza  L2.ancillary(i).cloud ...
                L2.ancillary(i).wind L2.ancillary(i).waves L2.ancillary(i).sal L2.ancillary(i).sst  L2.ancillary(i).binCount ...
                L2.Reflectance.Rrs(i,:) L2.Reflectance.Rrs_unc(i,:)
                ];
        end
        outFile.SB = strrep(outFile.SB,'ES_L2.sb','L2.sb'); %contingency
        outFile.SB = strrep(outFile.SB,'L2.sb','RRS_L2.sb');
    elseif strcmpi(dataTypes{dT},'Es')
        SBfields = 'station,date,time,lat,lon,RelAz,SZA,cloud,wind,waveht,sal,sst,bincount';
        SBunits = 'none,yyyyMMdd,hh:mm:ss,degrees,degrees,degrees,degrees,%,m/s,m,psu,degreesC,none';
        %% Es
        wavelength = L2.wavelength;
        for i=1:length(wavelength)
            SBfields = [SBfields  sprintf(',Es%.1f',wavelength(i))];
            SBunits = [SBunits ',uw/cm^2/nm'];
        end
        for i=1:length(wavelength)
            SBfields = [SBfields  sprintf(',Es%.1f_sd',wavelength(i))];
            SBunits = [SBunits ',uw/cm^2/nm'];
        end

        rows = length([L2.ancillary.datetime]);
        cols = ancCols + 2*length(wavelength);
        dataBlock = nan(rows,cols);
        for i=1:rows
            %% NOTE: Conversion from IMO L1B W/m^2/nm to uW/cm^2/nm
            dataBlock(i,:) = [
                L2.ancillary(i).station year(L2.ancillary(i).datetime) month(L2.ancillary(i).datetime) day(L2.ancillary(i).datetime) ...
                hour(L2.ancillary(i).datetime) minute(L2.ancillary(i).datetime) round(second(L2.ancillary(i).datetime)) ...
                L2.ancillary(i).lat L2.ancillary(i).lon L2.ancillary(i).relAz  L2.ancillary(i).sza  L2.ancillary(i).cloud ...
                L2.ancillary(i).wind L2.ancillary(i).waves L2.ancillary(i).sal L2.ancillary(i).sst  L2.ancillary(i).binCount ...
                L2.Irradiance.Es(i,:)/100 L2.Irradiance.Es_sd(i,:)/100
                ];
        end
        outFile.SB = strrep(outFile.SB,'RRS_L2.sb','L2.sb'); %contingency
        outFile.SB = strrep(outFile.SB,'L2.sb','ES_L2.sb');
    end

    dataBlock(isnan(dataBlock)) = nullChar;

    fidOut = fopen(outFile.SB,'w');
    fprintf(fidOut,'/begin_header\n');
    % Header
    for i=1:length(hFields)
        if ~strcmp(hFields{i}, 'fields') && ~strcmp(hFields{i}, 'units')
            if strcmp(hFields{i}, 'comments')
                for n=1:size(header.comments,1)
                    fprintf(fidOut,'!%s\n',header.comments(n,:));
                end
            else
                if ischar(header.(hFields{i}))
                    if contains(hFields{i},'time')
                        fprintf(fidOut,'/%s=%s[GMT]\n',hFields{i},header.(hFields{i}));
                    else
                        fprintf(fidOut,'/%s=%s\n',hFields{i},header.(hFields{i}));
                    end
                else
                    if contains(hFields{i},'latitude') || contains(hFields{i},'longitude')
                        fprintf(fidOut,'/%s=%.4f[DEG]\n',hFields{i},header.(hFields{i}));

                    else
                        fprintf(fidOut,'/%s=%.1f\n',hFields{i},header.(hFields{i}));
                    end
                end
            end
        end
    end
    fprintf(fidOut,'/instrument_manufacturer=IMO\n');
    fprintf(fidOut,'/instrument_model=DALEC\n');
    fprintf(fidOut,'/calibration_date=%s\n',sb.calibration_date);
    fprintf(fidOut,'/fields=%s\n',SBfields);
    fprintf(fidOut,'/units=%s\n',SBunits);
    fprintf(fidOut,'/end_header\n');

    % Body
    formatLine = '%.1f,%04d%02d%02d,%02d:%02d:%02d,%.4f,%.4f,%.1f,%.1f,%.1f,%.1f,%.1f,%.1f,%.1f,%d';
    for i=1:length(wavelength)*2
        formatLine = [formatLine ',%.6f'];
    end
    formatLine = [formatLine '\n'];

    for i=1:size(dataBlock,1)
        fprintf(fidOut, formatLine, dataBlock(i,:));
    end
    fclose(fidOut);
end