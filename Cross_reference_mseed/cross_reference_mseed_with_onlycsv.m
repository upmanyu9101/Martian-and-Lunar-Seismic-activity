function cross_reference_mseed_with_csv()
    % Select miniSEED files
    [filenames, pathname] = uigetfile('*.mseed', 'Select one or more miniSEED files', 'MultiSelect', 'on');
    
    if isequal(filenames, 0)
        disp('No files selected. Exiting...');
        return;
    end
    
    if ischar(filenames)
        mseedFiles = {fullfile(pathname, filenames)};
    else
        mseedFiles = fullfile(pathname, filenames);
    end

    % Select CSV file
    [csvFileName, csvPath] = uigetfile('*.csv', 'Select the CSV file for cross-referencing');
    if isequal(csvFileName, 0)
        disp('No CSV file selected. Exiting...');
        return;
    end
    
    csvFilePath = fullfile(csvPath, csvFileName);
    csvData = read_csv_data(csvFilePath); % Read CSV data
    
    % Loop through each miniSEED file
    for i = 1:length(mseedFiles)
        file = mseedFiles{i};
        fprintf('\nAnalyzing File: %s\n', file);

        % Read the miniSEED data
        [time, amplitude, metadata] = read_mseed_data(file);
        
        % Convert serial time to datetime format
        time = datetime(time, 'ConvertFrom', 'datenum', 'TimeZone', 'UTC');
        
        % Normalize amplitude for easier comparison
        amplitudeNorm = amplitude / max(abs(amplitude));
        
        % Cross-reference data at 5 points (start, end, and 3 in between)
        checkPoints = [1, round(length(time) * 0.25), round(length(time) * 0.5), ...
                       round(length(time) * 0.75), length(time)];
        
        matchSuccess = true;
        for k = 1:length(checkPoints)
            seismicOutput.time = time(checkPoints(k));
            seismicOutput.peakVelocity = max(abs(amplitudeNorm(checkPoints(k))));
            
            if ~cross_reference_with_csv(seismicOutput, csvData)
                matchSuccess = false;
                break;
            end
        end
        
        if matchSuccess
            fprintf('Data matches at all 5 points.\n');
        else
            fprintf('No match found or mismatch in CSV data.\n');
        end
    end
end

% Helper function to read CSV data
function csvData = read_csv_data(csvFilePath)
    csvData = readtable(csvFilePath);
end

% Helper function for cross-referencing with CSV data
function matchFound = cross_reference_with_csv(seismicOutput, csvData)
    matchFound = false;
    
    % Loop through CSV data to find matches
    for i = 1:height(csvData)
        csvTime = csvData.('time_abs__Y__m__dT_H__M__S__f_')(i);  % Replace with actual time column in CSV
        csvVelocity = csvData.('velocity_m_s_')(i);  % Replace with actual velocity column in CSV

        % Check if time and velocity match
        if abs(csvVelocity - seismicOutput.peakVelocity) < 0.01 && ...
           abs(datenum(seismicOutput.time) - datenum(csvTime)) < 1e-5
            matchFound = true;
            fprintf('Match found at time %s:\n', seismicOutput.time);
            disp(csvData(i, :));  % Display matching row
            break;
        end
    end
end

% Helper function to read miniSEED data
function [time, amplitude, metadata] = read_mseed_data(file)
    % This function reads miniSEED data from a file and returns time, amplitude, and metadata.

    % Read the miniSEED file using the rdmseed function (ensure the rdmseed.m is in your path)
    [seismicData, ~] = rdmseed(file);

    % Preallocate time and amplitude arrays
    numRecords = length(seismicData);
    time = [];
    amplitude = [];

    % Loop through each seismic data record and extract the time and amplitude
    for i = 1:numRecords
        time = [time; seismicData(i).t];  % Concatenate time data
        amplitude = [amplitude; seismicData(i).d];  % Concatenate amplitude data
    end

    % Extract metadata
    metadata.start_time = datestr(seismicData(1).RecordStartTimeMATLAB, 'yyyy-mm-ddTHH:MM:SS.FFF');
    metadata.sample_rate_hz = seismicData(1).SampleRate;
    metadata.station = seismicData(1).StationIdentifierCode;
    metadata.channel = seismicData(1).ChannelIdentifier;

end

