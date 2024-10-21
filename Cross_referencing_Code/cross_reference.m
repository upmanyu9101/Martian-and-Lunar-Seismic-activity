function cross_reference()
    % Allow the user to select both an mSEED file and a CSV file for comparison.
    % Prints matching time and velocity points between both datasets and prints the total matches at the end.

    % Select the MiniSEED file
    [mseedFile, mseedPath] = uigetfile('*.mseed', 'Select the MiniSEED file');
    if isequal(mseedFile, 0)
        disp('No MiniSEED file selected. Exiting...');
        return;
    end
    mseedFilePath = fullfile(mseedPath, mseedFile);

    % Read the MiniSEED data
    fprintf('Analyzing File: %s\n', mseedFilePath);
    [timeMseed, amplitudeMseed, metadata] = read_mseed_data(mseedFilePath);

    % Select the CSV file
    [csvFile, csvPath] = uigetfile('*.csv', 'Select the CSV file');
    if isequal(csvFile, 0)
        disp('No CSV file selected. Exiting...');
        return;
    end
    csvFilePath = fullfile(csvPath, csvFile);

    % Read the CSV data
    csvData = read_csv_data(csvFilePath);  % Read the CSV file into a table
    
    % Convert the time in MiniSEED data to datetime without timezone
    timeMseed = datetime(timeMseed, 'ConvertFrom', 'datenum');
    timeMseed.TimeZone = '';  % Remove time zone for comparison purposes

    % Define a tolerance for velocity matching
    velocityTolerance = 0;  % You can adjust this based on your data characteristics

    % Initialize match counter
    matchCount = 0;

    % Loop through each CSV time and match it with the closest MiniSEED time
    for i = 1:height(csvData)
        % Extract time and velocity from CSV
        csvTime = csvData.('time_abs__Y__m__dT_H__M__S__f_')(i);
         % Ensure this is the correct column name for time
        csvVelocity = csvData.('velocity_m_s_')(i);  % Ensure this is the correct column name for velocity

        % Find the closest MiniSEED time to the current CSV time
        [~, idx] = min(abs(timeMseed - csvTime));

        % Get the corresponding MiniSEED velocity
        mseedVelocity = amplitudeMseed(idx);

        % Check if the velocities match within the defined tolerance
        if abs(mseedVelocity - csvVelocity) == velocityTolerance
            % Print time and velocity from both sources in exponential notation
            fprintf('CSV Time: %s, CSV Velocity: %.5e\n', datestr(csvTime, 'yyyy-mm-dd HH:MM:SS.FFF'), csvVelocity);
            fprintf('MSEED Time: %s, MiniSEED Velocity: %.5e\n', datestr(timeMseed(idx), 'yyyy-mm-dd HH:MM:SS.FFF'), mseedVelocity);
            
            % Increment the match counter
            matchCount = matchCount + 1;
        end
    end

    % Print the total number of matches at the end
    fprintf('Total number of matching points: %d\n', matchCount);
end

% Helper function to read MiniSEED data
function [time, amplitude, metadata] = read_mseed_data(file)
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

    % Additional metadata (if available)
    if isfield(seismicData(1), 'Azimuth')
        metadata.azimuth_deg = seismicData(1).Azimuth;
    end
    if isfield(seismicData(1), 'Dip')
        metadata.dip_deg = seismicData(1).Dip;
    end
end

% Helper function to read the CSV data
function csvData = read_csv_data(csvFilePath)
    csvData = readtable(csvFilePath);  % Read the CSV file into a table
end
