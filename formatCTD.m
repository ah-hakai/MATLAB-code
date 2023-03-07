%% MATLAB function to format processed CTD data from the Hakai database
% Currently written for CTD data with Hakai data quality flags

% *However*, the function works with non-flagged CTD files as
% well, although this approach has not been evaluated as extensively. 

% Version 1.2. Last updated March 7, 2023.

% USAGE: formatCTD(fname, formatOut)

% 1. 'fname' is the input file name with extension, either an .xlsx file
        % or a .csv file, in native downloaded format from the Hakai database. 

% 2. 'formatOut' is either '1' for an .xlsx file, or '2' for a .csv file.

    % .xlsx file output - slowest option but contains all information in a
        % single file. Output file contains 4 worksheets as follows:
        
    % 1st worksheet contains 'AV' (acceptable value) flagged data, with
        % one cast per day, chosen based on maximum pressure obtained during
        % the CTD casts, and then based on the earliest time of collection,
        % if multiple casts reach equal max pressure on the same date. 
    % 2nd worksheet contains the additional 'AV' flagged
        % data, representing the extra casts collected on survey dates already
        % included in the first workseet. 
    % 3rd worksheet contains the 'SVC' (suspect value, confirm) flagged
        % data, which has been processed identically to the AV data, but still
        % contains flag columns. This allows the user to evaluate the individual 
        % sensor data in SVC casts for usage. Note that *any* SVC flag in
        % any data variable results in the entire cast from that date being
        % separated from the AV dataset, so that there is likely good data
        % in many of the SVC casts.
    % 4th worksheet contains the processing notes that describe the
        % steps performed on the CTD data from the Hakai database,
        % including information on the TEOS-10 derived variables, the last
        % date of QC flagging, where RBR and SBE oxygen sensor data has been
        % combined, etc.
        
    % .csv file output - fastest option but produces three separate .csv files
        % that are equivalent to the worksheets in the .xlsx output file above. 
        % Does not currently output a file containing the processing notes

% Current CTD flagging scheme (Jan. 24, 2022): 
% AV = acceptable value, 
% SVC = suspect value, confirm
% SVD = suspect value, discard

% * Other notes: 

% 1. This function is dependent on the raw variable names and their order in
% the downloaded files from the Hakai database matching the code below. If 
% they've changed, edit the text below accordingly.

% 2. This function is slowest when using a large .xlsx file as input, and/or
% outputting a large .xlsx file. To speed up the function, use a .csv file as
% input and/or choose a .csv file as the output


function formatCTD(fname, formatOut)


%% 1. Validate arguments and read in data files

if nargin ~= 2
    warning('Incorrect number of arguments. formatCTD() requires 2 input arguments') 
    return
end
  
[~, ~, ext] = fileparts(fname); % Confirm input file format

if ~any(strcmp(ext, {'.xlsx', '.csv'}))
    warning('Unrecognized input file format. Please use an .xlsx or a .csv file with full filename.')
    return
end
    
% Confirm output file format
switch formatOut % could also choose a default output format
    case 1
        extOut = '.xlsx';
    case 2 
        extOut = '.csv';
    otherwise
        warning('Unrecognized argument for output file format')
        return
end

switch ext % Load CTD data file
    
    case '.xlsx'
    
    % Data worksheet
    opts = detectImportOptions(fname, 'Sheet', 2); % 2nd worksheet holds data variables in .xlsx downloads
    
    % Ensure all flag columns are formatted as strings, even if they're currently unused
    idxFlags = contains(opts.VariableNames, 'Flag'); 
    opts = setvartype(opts, idxFlags, 'string'); % Set all flag columns as strings (unused are set as numeric/double in Hakai database)
    
    % Ensure all measurable parameters are formatted as numbers, even if they're currently unused
    measuredDataVars = {'Depth_m_', 'Conductivity_mS_cm_', 'Temperature_degC_', 'Pressure_dbar_',... 
        'PAR_umolM_2S_1_', 'FluorometryChlorophyll_ug_L_', 'Turbidity_FTU_', 'DissolvedO2_mL_L_', ...
        'SBRINKODissolvedO2_mL_L_', 'Salinity_PSU_', 'SpecificConductivity'};
    opts = setvartype(opts, measuredDataVars, 'double');
    
    % Identify desired variables
    selectedDataVars = {'CastPK', 'Cruise', 'Station', 'Latitude', 'Longitude', 'MeasurementTime'...
        'Depth_m_', 'DepthFlag', 'Conductivity_mS_cm_', 'ConductivityFlag', 'Temperature_degC_',...
        'TemperatureFlag', 'Pressure_dbar_', 'PressureFlag', 'PAR_umolM_2S_1_', 'PARFlag',...
        'FluorometryChlorophyll_ug_L_', 'FluorometryChlorophyllFlag', 'Turbidity_FTU_', 'TurbidityFlag',...
        'DissolvedO2_mL_L_', 'DissolvedO2_mL_L_Flag', 'SBRINKODissolvedO2_mL_L_', 'SBRINKODissolvedO2Flag',...
        'Salinity_PSU_', 'SalinityFlag', 'SpecificConductivity', 'SpecificConductivityFlag'};
    opts.SelectedVariableNames  = selectedDataVars; % Import selected variables only
    data = readtable(fname, opts, 'ReadVariableNames', 1, 'Sheet', 2);
    
    % Meta worksheet
    opts = detectImportOptions(fname, 'Sheet', 1);
    
    selectedMetaVars = {'CastPK', 'CTDModel', 'Cruise', 'Station', 'StartTime'};
    opts.SelectedVariableNames = selectedMetaVars; % Import selected variables only
    meta = readtable(fname, opts, 'ReadVariableNames', 1, 'Sheet', 1);
    
    % Find duplicated column names
    sharedVars = ismember(meta.Properties.VariableNames, data.Properties.VariableNames);

    % Join metadata and data by their CastPK values only
    sharedVarnames = meta.Properties.VariableNames(sharedVars);
    df = join(data, meta, 'Keys', 'CastPK', 'KeepOneCopy', sharedVarnames); 
    
    case '.csv'
    
    opts = detectImportOptions(fname);
    
    % Ensure all flag columns are formatted as strings, even if they're currently unused
    idxFlags = contains(opts.VariableNames, 'Flag'); 
    opts = setvartype(opts, idxFlags, 'string'); % Set all flag columns as strings (unused are set as double in Hakai database)
    
    % Ensure all measurable parameters are formatted as numbers, even if they're currently unused
    measuredDataVars = {'Depth_m_', 'Conductivity_mS_cm_', 'Temperature_degC_', 'Pressure_dbar_',... 
        'PAR_umolM_2S_1_', 'FluorometryChlorophyll_ug_L_', 'Turbidity_FTU_', 'DissolvedO2_mL_L_', ...
        'SBRINKODissolvedO2_mL_L_', 'Salinity_PSU_', 'SpecificConductivity'};
    opts = setvartype(opts, measuredDataVars, 'double');
    df = readtable(fname, opts);  % Metadata already merged with data when downloaded in .csv format
end


if ~exist('proNotes', 'var')
    proNotes = readtable('Processing_notes_for_flagged_CTD_data.xlsx', 'ReadVariableNames', 0); % General description of process, separate file
else
    warning('Processing notes file unavailable');
end


%% 2. Discard CTD scans with unusable depth measurements, and *clear contents of SVD-flagged values*

% Discard scans with unusable depth values
df = df(df.Depth_m_ ~= -9.9900e-29,:); % Remove all scans with native-sensor flagged depth measurements. 
df = df(~strcmp(df.DepthFlag, 'SVD'),:); % Remove entire scans with no usable depth

% Replace all other SVD-flagged data with NaN
idxFlags = contains(df.Properties.VariableNames, 'Flag'); % Identify flag columns to search for SVD flags
flagColnames = df.Properties.VariableNames(idxFlags); % List of flag column names

idxSVD = strcmp('SVD', df{:,idxFlags}); % Get 2D index of all SVD flags in *subset of flag columns*
[rows, cols] = find(idxSVD); % Column coordinates of SVD flags in *flag subset*

SVDcolNames = flagColnames(unique(cols)); % Identify names of columns containing SVD flags
idxSVDdf = find(ismember(df.Properties.VariableNames, SVDcolNames)); % Column column indices *in the data frame* that contain SVD flags

% Set as missing all values associated with remaining SVD flags
df{rows, idxSVDdf-1} = missing; % Depends on flag columns *always* following associated data


%% 3. Derive TEOS-10 variables

% Ensure all lat/long values are filled, to permit GSW Toolbox TEOS-10
% calculations
idxLat = find(~ismissing(df.Latitude), 1, 'first'); % first case of non-missing latitude
df.Latitude(isnan(df.Latitude)) = df.Latitude(idxLat); % Set missing values to first non-missing lat record
df.Longitude(isnan(df.Longitude)) = df.Longitude(idxLat); % Assuming longitude is also available if latitude is 

% Absolute Salinity
df.SA = gsw_SA_from_SP(df.Salinity_PSU_, df.Pressure_dbar_, df.Longitude, df.Latitude);

% Conservative Temperature
df.CT = gsw_CT_from_t(df.SA, df.Temperature_degC_, df.Pressure_dbar_);

% Density and Potential Density referenced to 0 m
df.rho = gsw_rho_CT_exact(df.SA, df.CT, df.Pressure_dbar_);
df.rho_pot0 = gsw_rho_CT_exact(df.SA, df.CT, 0); 

%  Convert Dissolved Oxygen to units of umol/kg

    % *** First replace impossible oxygen values with NaN
    % O2 < 0 occurs from sensor problems
    % O2 > 15 (arbitrary threshold) occurs when sensor cap is left on
    % during CTD deployment    
    
    % * For .csv files, must convert class to double
    if iscell(df.SBRINKODissolvedO2_mL_L_)
        df.SBRINKODissolvedO2_mL_L_ = str2double(df.SBRINKODissolvedO2_mL_L_);
    end
    
    %  Replace all oxygen values below 0 with NaN
    df.DissolvedO2_mL_L_(df.DissolvedO2_mL_L_ < 0) = NaN; % SBE data
    df.SBRINKODissolvedO2_mL_L_(df.SBRINKODissolvedO2_mL_L_ < 0) = NaN; % RINKO data
    
    %  Replace all oxygen values above 15 mL/L with NaN
    df.DissolvedO2_mL_L_(df.DissolvedO2_mL_L_ > 15) = NaN; % SBE data
    df.SBRINKODissolvedO2_mL_L_(df.SBRINKODissolvedO2_mL_L_ > 15) = NaN; % RINKO data
    
O2umol_L = df.DissolvedO2_mL_L_./0.022391;
O2_umol_kg = O2umol_L./(df.rho./1000); % *** Note: using actual density, not potential    
df.O2_umol_kg = O2_umol_kg;

% Apparent Oxygen Utilization, in umol/kg
df.O2sol = gsw_O2sol(df.SA, df.CT, df.Pressure_dbar_, df.Longitude, df.Latitude);
df.AOU_umol_kg = df.O2sol - O2_umol_kg;


%% 4. Handle times and dates

% Define Time Zone, to allow datetime() to account for daylight savings time and compute UTC time.
df.StartTime_local = datetime(df.MeasurementTime, 'InputFormat', 'yyyy-MM-dd HH:mm:ss.SSS', 'TimeZone', 'America/Vancouver', 'Format', 'HH:mm:ss');

% Identify local start time, and derive UTC start time
df.StartTime_UTC = datetime(df.StartTime_local, 'TimeZone', 'UTC');

% Derive separate date components, and DOY
df.Year = year(df.StartTime_local); % This approach strips off time, permits cast separation by date, i.e. groups casts on same day
df.Month = month(df.StartTime_local);
df.Day = day(df.StartTime_local);
df.Collection_date = datetime(df.Year, df.Month, df.Day);
df.DOY = day(df.StartTime_local, 'dayofyear');


%% 5. Select and clean up output variable names

% Clean up variable names, remove trailing underscore
df.Properties.VariableNames = strip(df.Properties.VariableNames, 'right', '_');

% Separately identify SBE and RINKO oxygen values
df.Properties.VariableNames{strcmp(df.Properties.VariableNames, 'DissolvedO2_mL_L')} = 'O2_mL_L'; % SBE oxygen sensor data
df.Properties.VariableNames{strcmp(df.Properties.VariableNames, 'SBRINKODissolvedO2_mL_L')} = 'RINKO_O2_mL_L'; % RINKO oxygen sensor data


%% 6. Separate AV and SVC flagged data
   
% 6.1 Find all rows containing SVC flags, for any variable
flagcols = contains(df.Properties.VariableNames, 'Flag');
idxSVC = strcmp('SVC', df{:, flagcols});
[svc_rows, ~] = find(idxSVC); % row indices of all non-zero elements (note multiple repeats because idx is 2D)
[av_rows, ~] = find(idxSVC==0);
svc_rows = unique(svc_rows);
av_rows = unique(av_rows);

% 6.2 Subset all AV and SVC-flagged data

% Selected variables, include TEOS-10 outputs
av_vars = {'CastPK', 'Station', 'Latitude', 'Longitude', 'CTDModel',...
    'Collection_date', 'Year', 'Month', 'DOY', 'StartTime_local', 'StartTime_UTC'...
    'Depth_m', 'Conductivity_mS_cm', 'Temperature_degC', 'Pressure_dbar',...
    'PAR_umolM_2S_1', 'FluorometryChlorophyll_ug_L', 'Turbidity_FTU',...
    'O2_mL_L', 'RINKO_O2_mL_L', 'Salinity_PSU', 'SA', 'CT', 'rho',...
    'rho_pot0', 'O2_umol_kg', 'AOU_umol_kg', };

% Different order for SVC data, to pair flags with corresponding variables
svc_vars = {'CastPK', 'Station', 'Latitude', 'Longitude', 'CTDModel',...
    'Collection_date', 'Year', 'Month', 'DOY', 'StartTime_local', 'StartTime_UTC'...
    'Depth_m', 'DepthFlag', 'Conductivity_mS_cm', 'ConductivityFlag', ...
    'Temperature_degC', 'TemperatureFlag', 'Pressure_dbar', 'PressureFlag',...
    'PAR_umolM_2S_1', 'PARFlag', 'FluorometryChlorophyll_ug_L',...
    'FluorometryChlorophyllFlag', 'Turbidity_FTU', 'TurbidityFlag',...
    'O2_mL_L', 'DissolvedO2_mL_L_Flag', 'RINKO_O2_mL_L', 'SBRINKODissolvedO2Flag',...
    'Salinity_PSU', 'SalinityFlag', 'SA', 'CT', 'rho',...
    'rho_pot0', 'O2_umol_kg', 'AOU_umol_kg', };

df_svc = df(svc_rows, [df.Properties.VariableNames(flagcols), svc_vars]);
df_av = df(av_rows,ismember(df.Properties.VariableNames, av_vars)); % AV-flagged data, no flags columns

% Also find date of last QC'd data, for time-stamp in output file notes
% * Use entire data frame, to identify last flag in either av or svc datasets
qcIdx = ~ismissing(df.DepthFlag); % Index of non-missing QC flags in DepthFlag column

if sum(qcIdx) == 0 % If there are no Depth Flags at all (e.g. Stn. DFO2 data, as of Oct. 8, 2021)...
    lastQCscan = find(~ismissing(df.ConductivityFlag), 1, 'last'); % ...use Conductivity flags to determine last QC date
else
    lastQCscan = find(qcIdx, 1, 'last'); % Otherwise, use Depth Flag to get index of last scan QC'd.
end

lastQCdate = datetime(df.StartTime_local(lastQCscan), 'format', 'yyyy-MM-dd'); % Date of last QC'd scan
str_QCdate = cellstr(['Last date of QC-flagged data: ', char(lastQCdate)]);


%% 7. Re-order columns and sort by time, depth, and then cast PK

df_av_reordered = df_av(:, av_vars); % Matches existing datasets; uses order above
df_svc_reordered = df_svc(:, svc_vars); % Pairs flags with data
    
% Sort rows. By Start_time (1st), then by Depth (2nd)
df_svc_sorted = sortrows(df_svc_reordered, {'StartTime_local', 'Depth_m', 'CastPK'});
df_av_sorted = sortrows(df_av_reordered, {'StartTime_local', 'Depth_m', 'CastPK'});


%% 8. Select one cast per day from the AV-flagged data
% Mutiple casts on the same day at the same station

G = findgroups(df_av_sorted.Collection_date, df_av_sorted.Station);
selected_casts = splitapply(@oneCastPerDay, df_av_sorted.CastPK, df_av_sorted.Pressure_dbar, df_av_sorted.StartTime_local, G);

% Select unique CTD casts for each day from original dataframe
selected_av_data = df_av_sorted(ismember(df_av_sorted.CastPK, selected_casts),:);
additional_av_data = df_av_sorted(~ismember(df_av_sorted.CastPK, selected_casts),:);


%% 9. Write both AV and SV data to output file(s)

str = [df.Station{1}, '_CTD_' string(datestr(min(df_av.Collection_date), 'yyyymmdd')), '_to_', string(datestr(max(df_av.Collection_date), 'yyyymmdd'))];

switch extOut
    case '.csv'
    
        % Writes to a set of .csv files (fastest) 
        fnameOutAV = [strcat(str{:}),'_AV.csv'];
        fnameOutAV_additional = [strcat(str{:}),'_additional_AV.csv'];
        fnameOutSVC = [strcat(str{:}), '_SVC.csv'];
        writetable(selected_av_data, fnameOutAV);
        writetable(additional_av_data, fnameOutAV_additional);
        writetable(df_svc_sorted, fnameOutSVC);

    case '.xlsx'
        % Writes to an .xlsx file
        fnameOut = [strcat(str{:}),'.xlsx'];
        writetable(selected_av_data, fnameOut, 'Sheet', 'Acceptable_Value (AV) data', 'AutoFitWidth', 0); % Write acceptable data
        writetable(additional_av_data, fnameOut, 'Sheet', 'Extra AV data', 'AutoFitWidth', 0); % Write additional data
        writetable(df_svc_sorted, fnameOut, 'Sheet', 'Suspect_data_to_confirm (SVC)', 'AutoFitWidth', 0); % Write suspect data

        % Append flagging date stamp and processing notes
        writecell(str_QCdate, fnameOut, 'Sheet', 'Notes', 'Range', 'A1:A1');
        
        if exist('proNotes', 'var')
            writetable(proNotes, fnameOut, 'Sheet', 'Notes', 'Range', 'A3:B40', 'WriteVariableNames', 0, 'AutoFitWidth', 0); % CTD data processing notes
        end
end

end


%% 10. Select one cast per day
% Based first on max pressure, then earliest time of day if there's two or
% more casts with equal max pressure

function [selected_casts] = oneCastPerDay(castPK, pressure, time)

    k = find(pressure==max(pressure)); % Identify deepest cast by selecting the deepest CTD scan
    selected_casts = castPK(k);  

    if height(selected_casts) > 1
        selected_casts = castPK(time == min(time(k)));
    end
end
