function rdmseed_multi()
    % This function categorizes lunar seismic events based on waveform characteristics.
    % The user selects one or more seismic data files (in miniSEED format).
    
    % Open file selection dialog for miniSEED files
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
    
    % Define parameters
    lowFreqThreshold = 1.5;      % Threshold for low-frequency classification (Hz)
    highFreqThreshold = 5;       % High-frequency threshold (Hz)
    veryHighFreqThreshold = 10;  % Very high-frequency threshold (Hz)
    pWaveVelocity = 5.0;         % Estimated P-wave velocity in km/s (lunar crust)
    sWaveVelocity = 2.9;         % Estimated S-wave velocity in km/s (lunar crust)
    densityMoon = 2500;          % Estimated density of the lunar crust (kg/m³)
    
    % Loop through each miniSEED file
    for i = 1:length(mseedFiles)
        file = mseedFiles{i};
        fprintf('\nAnalyzing File: %s\n', file);
        
        % Read the mSEED data and extract time and amplitude
        [time, amplitude, metadata] = read_mseed_data(file);
        
        % Preprocess amplitude data
        amplitudeNorm = amplitude / max(abs(amplitude));
        
        % Perform time-domain analysis
        eventDuration = (time(end) - time(1));  % Duration in datetime format
        fprintf('Event Duration: %s\n', eventDuration);
        
        % Perform frequency domain analysis (FFT and PSD)
        fs = metadata.sample_rate_hz;  % Sampling rate in Hz
        [Pxx, F] = periodogram(amplitudeNorm, [], [], fs);  % Power spectral density
        [~, maxIdx] = max(Pxx);
        dominantFreq = F(maxIdx);
        fprintf('Dominant Frequency: %.2f Hz\n', dominantFreq);
        
        % Classify the event based on dominant frequency
        if dominantFreq < lowFreqThreshold
            eventType = 'Low Frequency (LF) Event';
        elseif dominantFreq >= lowFreqThreshold && dominantFreq < highFreqThreshold
            eventType = 'High Frequency (HF) Event';
        elseif dominantFreq >= highFreqThreshold && dominantFreq < veryHighFreqThreshold
            eventType = 'Very High Frequency (VF) Event';
        else
            eventType = 'Super High Frequency (SF) Event';
        end
        fprintf('Event Type: %s\n', eventType);
        
        % Amplitude scaling and velocity conversion
        if isfield(metadata, 'scale_factor')
    amplitudeScaled = amplitude / metadata.scale_factor;  % Apply scaling if available
else
    amplitudeScaled = amplitude;  % Use original amplitude if no scaling factor is present
        end

        peakVelocity = max(abs(amplitudeScaled));
        fprintf('Peak Velocity: %.5e m/s\n', peakVelocity);
        
        % Bandpass filter to remove noise
        filterLowCutoff = 0.1;  % Lower cutoff frequency for bandpass filter (Hz)
        filterHighCutoff = min(10, fs/2);  % Upper cutoff frequency capped at Nyquist (fs/2)
        filteredAmplitude = bandpass(amplitudeScaled, [filterLowCutoff filterHighCutoff], fs);
        
        noise = amplitudeScaled - filteredAmplitude;  % Noise is the difference between original and filtered signal
        
        % Calculate event energy (simplified energy estimate)
        eventEnergy = sum(filteredAmplitude .^ 2) * (1/fs);
        fprintf('Estimated Energy: %.5e Joules\n', eventEnergy);
        
        % Seismic Moment and Magnitude Calculation
        distanceEstimate = 100;  % Assume a default distance of 100 km (this can be refined)
        seismicMoment = calc_seismic_moment(peakVelocity, distanceEstimate, densityMoon, pWaveVelocity);
        momentMagnitude = calc_moment_magnitude(seismicMoment);
        fprintf('Seismic Moment: %.5e N·m\n', seismicMoment);
        fprintf('Moment Magnitude: %.2f\n', momentMagnitude);
        
        % Signal-to-Noise Ratio (SNR) Calculation
        snrValue = calc_snr(filteredAmplitude, time);
        fprintf('Signal-to-Noise Ratio (SNR): %.2f\n', snrValue);
        
        % Coda Wave Analysis (Attenuation and Q-Factor)
        codaQFactor = analyze_coda_wave(time, filteredAmplitude, fs);
        fprintf('Coda Q-Factor: %.2f\n', codaQFactor);
        
        % Spectral Ratio Calculation (S-Wave to P-Wave)
        spectralRatio = calc_spectral_ratio(Pxx, F, lowFreqThreshold, highFreqThreshold);
        fprintf('Spectral Ratio (S-Wave/P-Wave): %.2f\n', spectralRatio);
       
         % P-wave and S-wave analysis (to estimate event distance)
[pWaveArrival, sWaveArrival] = detect_p_s_waves(time, filteredAmplitude);  % Custom P and S wave detection

if ~isempty(pWaveArrival) && ~isempty(sWaveArrival)
    % Calculate time difference between S-wave and P-wave arrivals
    deltaT = sWaveArrival - pWaveArrival;  % This gives you a 'duration' object

    % Convert the duration to seconds
    deltaT_seconds = seconds(deltaT);  % Convert 'deltaT' to seconds

    % Calculate the estimated distance
    estimatedDistance = deltaT_seconds * pWaveVelocity * sWaveVelocity / (pWaveVelocity - sWaveVelocity);

    % Convert pWaveArrival and sWaveArrival to string format for printing
    pWaveArrivalStr = datestr(pWaveArrival, 'yyyy-mm-dd HH:MM:SS.FFF');
    sWaveArrivalStr = datestr(sWaveArrival, 'yyyy-mm-dd HH:MM:SS.FFF');

    % Now print the result using the correct format
    fprintf('P-wave Arrival: %s, S-wave Arrival: %s\n', pWaveArrivalStr, sWaveArrivalStr);
    
else
    fprintf('P-wave or S-wave not clearly detected.\n');
end


        
        % Plot the filtered event (seismic signal) and the noise
        plot_filtered_event(time, filteredAmplitude, noise, F, Pxx, eventType);
    end
end

% -- Helper Functions --

function seismicMoment = calc_seismic_moment(peakVelocity, distance, density, velocity)
    % Calculate the seismic moment (M₀)
    seismicMoment = peakVelocity * distance * velocity^3 * density;
end

function momentMagnitude = calc_moment_magnitude(seismicMoment)
    % Calculate moment magnitude (Mₜ) from seismic moment (M₀)
    momentMagnitude = (2/3) * log10(seismicMoment) - 6;
end

function snrValue = calc_snr(amplitude, time)
    % Calculate Signal-to-Noise Ratio (SNR)
    noiseRegion = amplitude(1:round(0.1 * length(time)));  % Assume the first 10% is noise
    signalRegion = amplitude(round(0.1 * length(time)):end);  % The rest is the signal
    snrValue = 20 * log10(rms(signalRegion) / rms(noiseRegion));
end

function codaQFactor = analyze_coda_wave(time, amplitude, fs)
    % Coda wave analysis and Q-factor estimation
    codaStart = round(0.8 * length(time));  % Assume coda starts at 80% of the waveform
    codaWave = amplitude(codaStart:end);
    codaEnergy = sum(codaWave.^2);
    
    % Convert the coda duration to seconds
    codaDuration = length(codaWave) / fs;  % Duration in seconds
    
    % Q-factor estimate based on energy and duration
    codaQFactor = codaEnergy / codaDuration;  % Simplified Q-factor estimate
end

function spectralRatio = calc_spectral_ratio(Pxx, F, lowFreq, highFreq)
    % Calculate the spectral ratio between S-wave and P-wave regions
    pWaveEnergy = sum(Pxx(F < lowFreq));
    sWaveEnergy = sum(Pxx(F > lowFreq & F < highFreq));
    spectralRatio = sWaveEnergy / pWaveEnergy;
end

function [pWaveArrival, sWaveArrival] = detect_p_s_waves(time, amplitude)
    % Assume pWaveArrival is precomputed or detected in some way
    pWaveArrival = 1000;  % Example index for P-wave arrival (you may have a better method)

    % Reduce memory load by processing smaller segments of amplitude
    batch_size = 10000;  % Example batch size
    sWaveArrival = [];

    for start_idx = 1:batch_size:length(amplitude)
        end_idx = min(start_idx + batch_size - 1, length(amplitude));
        current_batch = amplitude(start_idx:end_idx);
        
        % Detect S-wave within the current batch
        max_amplitude = max(current_batch);
        sWave_idx = find(current_batch >= 0.1 * max_amplitude & (start_idx:end_idx) > pWaveArrival, 1, 'first');
        
        if ~isempty(sWave_idx)
            sWaveArrival = start_idx + sWave_idx - 1;
            break;
        end
    end

    if isempty(sWaveArrival)
        warning('S-wave arrival not detected.');
    end
end


function plot_filtered_event(time, filteredAmplitude, noise, F, Pxx, eventType)
    % Plot the filtered seismic event (without noise) and the noise separately
    
    figure;
    
    % Plot the filtered waveform (seismic event)
    subplot(3, 1, 1);
    plot(time, filteredAmplitude);
    xlabel('Time (UTC)');
    ylabel('Amplitude (filtered)');
    title(['Filtered Waveform - ', eventType]);
    
    % Plot the noise (removed from the original signal)
    subplot(3, 1, 2);
    plot(time, noise);
    xlabel('Time (UTC)');
    ylabel('Noise Amplitude');
    title('Noise (Removed from Signal)');
    
    % Plot the power spectral density
    subplot(3, 1, 3);
    plot(F, 10 * log10(Pxx));
    xlabel('Frequency (Hz)');
    ylabel('Power/Frequency (dB/Hz)');
    title('Power Spectral Density');
    
    grid on;
end

function [time, amplitude, metadata] = read_mseed_data(file)
    % This function reads miniSEED data from a file and returns time, amplitude, and metadata.
    
    % Read the miniSEED file using rdmseed function
    [seismicData, ~] = rdmseed(file);
    % Preallocate time and amplitude arrays
numRecords = length(seismicData);
time = []; % Alternatively: preallocate if the total size is known
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

    
    % Additional metadata (if available from seismicData or channelInfo)
    if isfield(seismicData, 'Azimuth')
        metadata.azimuth_deg = seismicData(1).Azimuth;
    end
    if isfield(seismicData, 'Dip')
        metadata.dip_deg = seismicData(1).Dip;
    end
    
    % You can expand this section to add more metadata fields as necessary
end
