function ngc_4e5e
    global isMarsData;

    % Ask user for input to toggle between Mars and lunar data
    choice = input('Analyze Mars data or Lunar data? Enter "m" for Mars, "l" for Lunar: ', 's');
    
    if strcmp(choice, 'm')
        isMarsData = true;  % Analyze Mars data
        disp('Analyzing Mars data.');
    elseif strcmp(choice, 'l')
        isMarsData = false;  % Analyze Lunar data
        disp('Analyzing Lunar data.');
    else,
        disp('Invalid input. Defaulting to Lunar data.');
        isMarsData = false;  % Default to Lunar data
    end

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
    
    % Loop through each miniSEED file
    for i = 1:length(mseedFiles)
        file = mseedFiles{i};
        fprintf('\nAnalyzing File: %s\n', file);
        
        % Read the mSEED data and extract time and amplitude
        [time, amplitude, metadata] = read_mseed_data(file);
        
        % Convert serial time to datetime format
        time = datetime(time, 'ConvertFrom', 'datenum', 'TimeZone', 'UTC');
        
        % Preprocess amplitude data
        amplitudeNorm = amplitude / max(abs(amplitude));

        % Detect the seismic event based on amplitude threshold
        noiseThreshold = 0.01;  % Adjust this value based on your data characteristics
        [eventStartIdx, eventEndIdx] = detect_event(amplitudeNorm, noiseThreshold);
        
        % Redefine the time and amplitude for the actual seismic event
        eventTime = time(eventStartIdx:eventEndIdx);
        eventAmplitude = amplitudeNorm(eventStartIdx:eventEndIdx);
        
        % Calculate event duration based on detected event
        eventDuration = eventTime(end) - eventTime(1);  % Calculate event duration
        fprintf('Event Duration (Waveform Duration): %s\n', eventDuration);


        % Perform frequency domain analysis (FFT and PSD)
        fs = metadata.sample_rate_hz;  % Sampling rate in Hz
        [Pxx, F] = periodogram(amplitudeNorm, [], [], fs);  % Power spectral density
        [~, maxIdx] = max(Pxx);
        dominantFreq = F(maxIdx);
        fprintf('Dominant Frequency: %.2f Hz\n', dominantFreq);
        
         % Classify the event based on dominant frequency
        eventType = classify_event(dominantFreq);
        fprintf('Event Type: %s\n', eventType);

        % Choose MarsF or lunar parameters based on user input
        if isMarsData
            % Parameters for Mars
            pWaveVelocity = 3.5;         % P-wave velocity for Martian crust (km/s)
            sWaveVelocity = 2.0;         % S-wave velocity for Martian crust (km/s)
            density = 2900;              % Martian crust density (kg/m³)
            fprintf('Using Mars parameters.\n');
        else
            % Parameters for Moon
            pWaveVelocity = 5.0;         % P-wave velocity for lunar crust (km/s)
            sWaveVelocity = 2.9;         % S-wave velocity for lunar crust (km/s)
            density = 2500;              % Lunar crust density (kg/m³)
            fprintf('Using lunar parameters.\n');
        end
        
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
        filterHighCutoff = dynamic_upper_cutoff(fs);  % Adjusted upper cutoff based on the data
        filteredAmplitude = bandpass(amplitudeScaled, [filterLowCutoff filterHighCutoff], fs);
        
        noise = amplitudeScaled - filteredAmplitude;  % Noise is the difference between original and filtered signal
        
        % Calculate event energy (simplified energy estimate)
        eventEnergy = sum(filteredAmplitude .^ 2) * (1/fs);
        fprintf('Estimated Energy: %.5e Joules\n', eventEnergy);
        
        % Seismic Moment and Magnitude Calculation
        distanceEstimate = 100;  % Assume a default distance of 100 km (this can be refined)
        seismicMoment = calc_seismic_moment(peakVelocity, distanceEstimate, density, pWaveVelocity);
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
        spectralRatio = calc_spectral_ratio(Pxx, F, 1.5, 5);  % Example frequency ranges
        fprintf('Spectral Ratio (S-Wave/P-Wave): %.2f\n', spectralRatio);

        % P-wave and S-wave analysis (using CWT)
        [pWaveIdx, sWaveIdx] = detect_p_s_waves_cwt(time, amplitudeNorm, fs);

        % Check if both P-wave and S-wave indices were detected
        if ~isempty(pWaveIdx) && ~isempty(sWaveIdx)
            pWaveArrival = time(pWaveIdx);
            sWaveArrival = time(sWaveIdx);
            fprintf('P-wave Arrival: %s, S-wave Arrival: %s\n', datestr(pWaveArrival, 'yyyy-mm-dd HH:MM:SS.FFF'), datestr(sWaveArrival, 'yyyy-mm-dd HH:MM:SS.FFF'));
            
            % Calculate estimated distance
            deltaT = sWaveArrival - pWaveArrival;
            deltaT_seconds = seconds(deltaT);
            estimatedDistance = deltaT_seconds * pWaveVelocity * sWaveVelocity / (pWaveVelocity - sWaveVelocity);
            fprintf('Estimated Distance: %.2f km\n', estimatedDistance);
        else
            fprintf('P-wave or S-wave not clearly detected.\n');
        end
        
        % Plot the filtered event (seismic signal) and the noise
        plot_filtered_event(time, filteredAmplitude, noise, F, Pxx, 'Seismic Event', pWaveArrival, sWaveArrival);
    end
end


% Callback function to toggle between Mars and lunar data
function keyPressCallback(~, event)
    global isMarsData;
    if strcmp(event.Key, 'm')  % Press "m" to toggle between Mars and lunar data
        isMarsData = ~isMarsData;  % Toggle the Mars/Lunar state
        if isMarsData
            disp('Now analyzing Mars data.');
        else
            disp('Now analyzing lunar data.');
        end
    end
end

% Helper to Read MiniSEED Data
function [time, amplitude, metadata] = read_mseed_data(file)
    [seismicData, ~] = rdmseed(file);  % Read miniSEED file using rdmseed
    numRecords = length(seismicData);
    time = [];
    amplitude = [];
    for i = 1:numRecords
        time = [time; seismicData(i).t];  % Concatenate time data
        amplitude = [amplitude; seismicData(i).d];  % Concatenate amplitude data
    end
    metadata.start_time = datestr(seismicData(1).RecordStartTimeMATLAB, 'yyyy-mm-ddTHH:MM:SS.FFF');
    metadata.sample_rate_hz = seismicData(1).SampleRate;
    metadata.station = seismicData(1).StationIdentifierCode;
    metadata.channel = seismicData(1).ChannelIdentifier;

    if isfield(seismicData(1), 'Azimuth')
        metadata.azimuth_deg = seismicData(1).Azimuth;
    end
    if isfield(seismicData(1), 'Dip')
        metadata.dip_deg = seismicData(1).Dip;
    end
end

% Event Detection Helper
function [eventStartIdx, eventEndIdx] = detect_event(signal, threshold)
    eventIndices = find(abs(signal) > threshold);
    if isempty(eventIndices)
        eventStartIdx = 1;  % Default to start of the data
        eventEndIdx = length(signal);  % Default to end of the data
    else
        eventStartIdx = eventIndices(1);  % First point where amplitude exceeds threshold
        eventEndIdx = eventIndices(end);  % Last point where amplitude exceeds threshold
    end
end

% P-wave and S-wave Detection Helper using CWT
function [pWaveIdx, sWaveIdx] = detect_p_s_waves_cwt(time, amplitude, fs)
    % Perform Continuous Wavelet Transform (CWT) on the amplitude signal
    [wt, f] = cwt(amplitude, 'amor', fs);

    % Define frequency bands for P-wave and S-wave based on expected ranges
    pWaveFreqRange = [1, 5];  % Example range for P-wave frequency
    sWaveFreqRange = [0.5, 3];  % Example range for S-wave frequency

    % Sum across frequency bands to emphasize specific ranges
    pWaveEnergy = sum(abs(wt(f >= pWaveFreqRange(1) & f <= pWaveFreqRange(2), :)), 1);  % Summing across the P-wave frequency range
    sWaveEnergy = sum(abs(wt(f >= sWaveFreqRange(1) & f <= sWaveFreqRange(2), :)), 1);  % Summing across the S-wave frequency range

    % Detect peaks in the P-wave and S-wave energy
    [~, pWaveIdx] = max(pWaveEnergy);  % P-wave arrival is the peak in the P-wave energy
    [~, sWaveIdx] = max(sWaveEnergy(pWaveIdx+1:end));  % S-wave arrival is the peak after the P-wave
    sWaveIdx = sWaveIdx + pWaveIdx;  % Adjust S-wave index relative to the full signal
end

% Seismic Moment and Magnitude Calculation Helper
function seismicMoment = calc_seismic_moment(peakVelocity, distance, density, velocity)
    seismicMoment = peakVelocity * distance * velocity^3 * density;
end

% Moment Magnitude Calculation Helper
function momentMagnitude = calc_moment_magnitude(seismicMoment)
    momentMagnitude = (2/3) * log10(seismicMoment) - 6;
end

% Signal-to-Noise Ratio (SNR) Calculation Helper
function snrValue = calc_snr(amplitude, time)
    noiseRegion = amplitude(1:round(0.1 * length(time)));  % Assume the first 10% is noise
    signalRegion = amplitude(round(0.1 * length(time)):end);  % The rest is the signal
    snrValue = 20 * log10(rms(signalRegion) / rms(noiseRegion));
end

% Coda Wave Analysis (Attenuation and Q-Factor) Helper
function codaQFactor = analyze_coda_wave(time, amplitude, fs)
    codaStart = detect_coda_start(amplitude);
    if isempty(codaStart)
        codaQFactor = 0;
        return;
    end
    
    codaWave = amplitude(codaStart:end);
    codaEnergy = sum(codaWave.^2);
    codaDuration = length(codaWave) / fs;  % Duration in seconds

    if codaEnergy < 1e-10 || codaDuration == 0
        codaQFactor = 0;
    else
        codaQFactor = codaEnergy / codaDuration;
    end
end

% Coda Wave Start Detection Helper
function codaStart = detect_coda_start(amplitude)
    threshold = 0.01 * max(amplitude);  % Example threshold as 1% of max amplitude
    codaStart = find(amplitude < threshold, 1, 'first');
    if isempty(codaStart)
        codaStart = [];
    end
end

% Spectral Ratio Calculation Helper
function spectralRatio = calc_spectral_ratio(Pxx, F, lowFreq, highFreq)
    pWaveEnergy = sum(Pxx(F < lowFreq));
    sWaveEnergy = sum(Pxx(F > lowFreq & F < highFreq));
    spectralRatio = sWaveEnergy / pWaveEnergy;
end

% Dynamic Upper Cutoff Frequency Helper
function filterHighCutoff = dynamic_upper_cutoff(fs)
    nyquistFreq = fs / 2;
    filterHighCutoff = min(nyquistFreq * 0.8, 10);  % Max upper cutoff at 10 Hz or less
end

% Filtered Event Plotting Helper
function plot_filtered_event(time, filteredAmplitude, noise, F, Pxx, eventType, pWaveArrival, sWaveArrival)
    % Plot the filtered waveform
    subplot(3, 1, 1);
    plot(time, filteredAmplitude);
    hold on;
    if ~isempty(pWaveArrival)
        xline(pWaveArrival, '--r', 'P-wave Arrival', 'LabelVerticalAlignment', 'bottom', 'LabelHorizontalAlignment', 'center');
    end
    if ~isempty(sWaveArrival)
        xline(sWaveArrival, '--b', 'S-wave Arrival', 'LabelVerticalAlignment', 'bottom', 'LabelHorizontalAlignment', 'center');
    end
    xlabel('Time (UTC)');
    ylabel('Amplitude (filtered)');
    title(['Filtered Waveform - ', eventType]);
    grid on;
    hold off;

    % Plot the noise (removed from the original signal)
    subplot(3, 1, 2);
    plot(time, noise);
    xlabel('Time (UTC)');
    ylabel('Noise Amplitude');
    title('Noise (Removed from Signal)');
    grid on;

    % Plot the power spectral density
    subplot(3, 1, 3);
    plot(F, 10 * log10(Pxx));
    xlabel('Frequency (Hz)');
    ylabel('Power/Frequency (dB/Hz)');
    title('Power Spectral Density');
    grid on;
end
% Event Classification Helper
function eventType = classify_event(dominantFreq)
    % Define the frequency thresholds for classification
    lfThreshold = 0.5;  % Low Frequency (LF) threshold
    hfThreshold = 5;    % High Frequency (HF) threshold
    vfThreshold = 10;   % Very High Frequency (VF) threshold
    sfThreshold = 20;   % Super High Frequency (SF) threshold
    resonanceThreshold = 30;  % Example resonance frequency threshold

    % Classify based on frequency ranges
    if dominantFreq < lfThreshold
        eventType = 'Low Frequency (LF) Event';
    elseif dominantFreq >= lfThreshold && dominantFreq < hfThreshold
        eventType = 'High Frequency (HF) Event';
    elseif dominantFreq >= hfThreshold && dominantFreq < vfThreshold
        eventType = 'Very High Frequency (VF) Event';
    elseif dominantFreq >= vfThreshold && dominantFreq < sfThreshold
        eventType = 'Super High Frequency (SF) Event';
    else
        eventType = 'Resonance Event';
    end
end

