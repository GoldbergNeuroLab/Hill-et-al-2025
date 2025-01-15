%% Analysis of intrinsic properties and EAP firing from current steps
% Script Author: Sophie Hill, adapted from Sophie Liebergall
% Updated: 2 February 2024
% Inputs:
%   single .abf file with trace from current steps protocol in IC (5 pA
%   steps startging from -100 pA)
% Outputs:
%   Intrinsic properties
%   EAP properties
% Dependencies:
%   - requires abfload.m and extrema.m to be in MATLAB search path

%% Import File

% Open a dialog box to select the file of interest from all .abf files in
% the folder
% [fileName,path] = uigetfile('*.abf', 'abf files');
% filePath = [path fileName];
% make this the default path
% cd(path)

[data,si,~]=abfload(filePath);
% data = data frame configured as <data pts per sweep> by <number of chans> 
%   by <number of sweeps>.
% si = sampling interval in microseconds (us)
%       Note: sampling rate (Hz) = 1e6 / si (us)
% h = information on file (selected header parameters)

% Extract dimensions of the full data matrix
[i, num_chan, num_sweeps]=size(data); 
% i = number of samples in each sweep (i.e. recorded data points)
% num_chan = number of channels
% num_sweeps = number of sweeps in .abf file

% Specify channel of interest (assuming first channel is membrane voltage
% and second is current)
% User may need to change this value depending on channel of interest
v_chan = 1;
i_chan = 2;
    
% Create vector with all indices in file
indices = 1:i;
clear i

% Create a vector with all indices converted to time points
time = indices ./ (1e6/si);

% Find on and off times for current steps
I_diff = diff(data(:,i_chan,1)); % extract current signal (channel 2) for first sweep, then take the first derivative at each point
I_on = find(I_diff == min(I_diff),1); % get index when current goes on (minimum of dI/dt)
I_off = find(I_diff == max(I_diff),1); % get index when current goes off (maximum of dI/dt)
% Note: this assumes that this time is the same for each sweep and that the
% first sweep injects negative current
clear I_diff indices path h
%% Create Table of Spikes for Each Sweep

% Set threshold for defining an action potential (in mV)
dv_dt_threshold = 10;
v_threshold = -10;
double_counting = 0.0008 * (1e6/si);

spike_table = zeros(num_sweeps,3); % create empty table to store spikes for each sweep
% Formatted as [Sweep] [Absolute Current (pA)] [Difference from Holding (pA)] [Number of Spikes] 

spike_i = {}; % create empty cell array to store indices of spikes for each sweep (to calc ISIs)
dt = (si*1000) / 1e6; % units ms
spike_x = {};
spike_y = {};

for sweep = 1:num_sweeps % loop through each sweep in the .abf file
    v = data(I_on:I_off+100,v_chan,sweep); % extract voltage data (channel 1) for current sweep
    dv = gradient(v); % dV, units mV
    dv_dt = dv ./ dt; % dV/dt, units mV/ms

    k_v = find(v > v_threshold);
    k_dv = find(dv_dt > dv_dt_threshold);
    k = intersect(k_v, k_dv);

    spike_xs = [];
    spike_ys = [];

    if isempty(k)
        spike_is = [];
    else
        spike_is = [k(1)];
        spike_xs = [spike_xs, (k(1) + I_on) / (1e6/si)];
        spike_ys = [spike_ys, v(k(1))];
        for w = 2:length(k)
            if (k(w) - k(w-1)) > 300
                spike_is = [spike_is, k(w)];
                spike_xs = [spike_xs, (k(w) + I_on) / (1e6/si)];
                spike_ys = [spike_ys, v(k(w))];
            end
        end
    end    
    spike_i{sweep} = spike_is; % get indices of local maxima above threshold (i.e. spikes)
    spike_x{sweep} = spike_xs;
    spike_y{sweep} = spike_ys;

    I = data(int16((I_on+I_off)/2),i_chan,sweep); % get current value for sweep
    I_hold = data(int16(I_on-100),i_chan,sweep);

    spike_table(sweep,1) = sweep; % record sweep #
    spike_table(sweep,2) = I; % record absolute current for sweep
    spike_table(sweep,3) = I - I_hold; % record current difference from hold for sweep
    spike_table(sweep,4) = size(spike_i{sweep},2); % record number of spikes for a sweep
end
clear v dv dv_dt k w spike_is sweep I I_hold k_dv k_v dt

spike_sweeps = spike_table(spike_table(:,4) > 0);
rheo_sweep = spike_sweeps(1);

%% Passive membrane properties

% Resting membrane potential
abs_current = (abs(spike_table(:,3)));
idx_abs_current = find(abs_current == min(abs_current));
Vm = mean(data(I_on:I_off,v_chan,idx_abs_current));
clear abs_current idx_abs_current

% Input Resistance, calculated from most negative current injection
neg_sweeps = spike_table(spike_table(:,3) < 0);
sort_neg_sweep = sort(neg_sweeps);
most_neg_sweep = sort_neg_sweep(1);
V1 = mean(data(1:I_on,v_chan,most_neg_sweep)); % Calcualte mean voltage before current step
V2 = mean(data(((I_on+I_off)./2):I_off,v_chan,most_neg_sweep)); % Calclate mean voltage in second half of current step
I1 = mean(data(1:I_on,i_chan,most_neg_sweep)); % Calculate mean current before current step
I2 = mean(data(((I_on+I_off)./2):I_off,i_chan,most_neg_sweep)); % Calclate mean current in second half of current step
Rm = (V2 - V1)./(I2 - I1).*1000; % Calculate Rm in MOhms
clear V1 V2 I1 I2 neg_sweeps sort_neg_sweep 

% Membrane time constant, calculated from exponential fit of 120 ms 
fit_dur = 20; % units ms
fit_data = data(I_on:(I_on+(fit_dur.*1e3./si)),v_chan,most_neg_sweep);
fit_data_time = transpose((1:fit_dur*1e3./si+1)./1e3.*si); % Create time bin for calculating membrane time constant from first sweep
ft = fittype('a*exp(-b*t) + c','indep','t'); % Specify fit type as single expontential with c term
start_pt = [6,0.1,-65]; % Set start point for fit coefficients (from literature) [6,0.1,-65]
f = fit(fit_data_time,fit_data,ft,'StartPoint',start_pt); % Fit curve
% figure;
% plot(f,fit_data_time,fit_data)
membrane_time_constant = 1/f.b;
clear fit_dur ft start_pt ft

% Sag Percentage
sag_min = min(data(I_on:I_off,v_chan,most_neg_sweep)); % Get min voltage during negative current inj.
steady_index = int16((I_off - I_on)*0.75 + I_on); % Calculate steady state as voltage 3/4 of way through current step
sag_steady = data(steady_index,v_chan,most_neg_sweep);
sag = 100 * (sag_steady/sag_min); % Units percent
clear sag_min steady_index sag_steady most_neg_sweep

%% Properties of tonic firing

% Max Steady State Firing Frequency
[max_spikes, ~] = max(transpose(spike_table(:,4))); % Select the sweep with the greatest number of spikes
max_SSFF = max_spikes./((I_off - I_on)./1e6.*si); % in Hz
clear max_spikes

% Max Instantaneous Firing Frequency
ISIs = {}; % Create cell array with inter-spike intervals (ISIs)
for sweep = 1:num_sweeps % loop through each sweep
    ISI_sweep = []; % create vector to store ISIs for each sweep
    if spike_table(sweep,4) > 1
        for int = 2:spike_table(sweep,4) % loop through each ISI
            % Calculate ISI and append to ISI_sweep vector
            ISI_sweep(int-1) = [(spike_i{sweep}(int) - spike_i{sweep}(int-1))];
        end
    ISIs{sweep} = ISI_sweep; % Add ISIs for a sweep to ISIs cell array
    end
end
sorted_ISIs = sort(cell2mat(ISIs)); 
max_IFF = 1 ./ (sorted_ISIs(1) ./ 1e6 .* si);
clear ISI_sweep sweep int sorted_ISIs

% Define tonic sweep
ton_sweeps = spike_table(spike_table(:,4) > 40);
if ~isempty(ton_sweeps)
    tonic_sweep = ton_sweeps(1);
else 
    ton_sweeps = spike_table(spike_table(:,4) == max(spike_table(:,4)));
    tonic_sweep = ton_sweeps(1);
end
clear ton_sweeps

% Spike frequency adaptation
sfa_1_v_2 = (ISIs{tonic_sweep}(1) / ISIs{tonic_sweep}(2));
if spike_table(tonic_sweep, 4) >= 11
    sfa_1_v_10 = (ISIs{tonic_sweep}(1) / ISIs{tonic_sweep}(10));
else sfa_1_v_10 = NaN;
end 
sfa_1_v_n = (ISIs{tonic_sweep}(1) / ISIs{tonic_sweep}(end));

% Properties from first spike at rheobase

% Rheobase
spike_sweeps = spike_table(spike_table(:,4) > 0);
rheo_sweep = spike_sweeps(1); % Get sweep number of rheobase (rheo_sweep)
rheobase = spike_table(rheo_sweep,3); % Get current injection (pA) during rheo_sweep sweep (difference from hold)
clear spike_sweeps

% Spikes at rheobase
spikes_at_rheobase = spike_table(rheo_sweep,4);

% AP threshold at rheobase
endtime = round(0.003*(1e6./si));
if length(spike_i{rheo_sweep}) == 1
    otherside = 0.1*(1e6./si);
else
    otherside = spike_i{rheo_sweep}(2) - spike_i{rheo_sweep}(1) - endtime;
end
AP1_v = data(I_on:(I_on+spike_i{rheo_sweep}(1)+otherside),v_chan,rheo_sweep);
AP1_v_d = gradient(AP1_v)./(si*1e-3); % get 1st derivative of voltage data (mV/ms)
AP_thresh_i = find(AP1_v_d(100:end) >= 10);
AP_thresh_i = AP_thresh_i(3)+99; % Get index where dV/dt >/ 10 mV/ms
AP_thresh_rheobase = AP1_v(AP_thresh_i);
clear otherside

% Latency to first spike at rheobase
latency_to_AP1_rheobase = AP_thresh_i./1e3*si;

% AP Peak
[AP_peak, AP_peaki] = max(AP1_v);

% AP Rise Time, calculated as difference between time at AP peak and time at AP threshold
AP_rise_time = (AP_peaki - AP_thresh_i)./1e3.*si; %in ms

% Max upstroke velocity, calculated as max of 1st derivative between AP threshold and peak
max_rise_slope = max(AP1_v_d(AP_thresh_i:AP_peaki));
clear AP_thresh_i 

% AP Amplitude
AP_amplitude = AP_peak - AP_thresh_rheobase;

% AP Halfwidth, calculated using AP threshold and AP peak. 
AP_half_up = find(AP1_v >= (AP_amplitude*0.5+AP_thresh_rheobase));
AP_half_up = AP_half_up(1); % Find index of 1/2 AP amplitude on upstroke
AP_half_down = find(AP1_v(AP_peaki:end) <= (AP_amplitude*0.5+AP_thresh_rheobase));
AP_half_down = AP_half_down(1) + AP_peaki - 1; % Find index of 1/2 AP amplitude on downstroke
AP_half_width = (AP_half_down - AP_half_up)./1e3*si;
clear AP_half_up AP_half_down

% After Hyperpolarization (AHP) Amplitude, calculated as difference between AP threshold and AHP trough V
[AHP_min, AHP_min_time] = min(AP1_v(AP_peaki:end));
AHP_amplitude = AP_thresh_rheobase - AHP_min;
clear AHP_min 

% Max downstroke velocity
max_decay_slope = min(AP1_v_d(AP_peaki:(AP_peaki + AHP_min_time-1)));
clear AP1_v_d AP_peaki AHP_min_time

%% Identify large EAPs (dV/dt > 10 & v > -10)
eap_table = zeros(num_sweeps,3); % create empty table to store spikes for each sweep
% Formatted as:
% Col 1 = sweep number
% Col 2 = absolute current injection (pA)
% Col 3 = difference from holding current (pA)
% Col 4 = number of EAPs

eap_i = {}; % create empty cell array to store indices of spikes for each sweep (to calc ISIs)
dt = (si*1000) / 1e6; % units ms
large_x = {}; % x indices for plotting EAPs
large_y = {}; % y indices for plotting EAPs
eap_v_threshold = -20;

for sweep = 1:num_sweeps % loop through each sweep in the .abf file
    v = data(I_off+(0.008*(1e6/si)):end,v_chan,sweep); % extract voltage data (channel 1) after current step
    dv = gradient(v); % dV, units mV
    dv_dt = dv ./ dt; % dV/dt, units mV/ms
    
    k_v = find(v > eap_v_threshold);
    k_dv = find(dv_dt > dv_dt_threshold);
    k = intersect(k_v, k_dv);

    large_yw = [];
    large_xw = [];
    if isempty(k)
        eap_is = [];
    else
        eap_is = [k(1)];
        large_yw = [large_yw, data(k(1) + I_off+(0.008*(1e6/si)), v_chan, sweep)];
        large_xw = [large_xw, (k(1) + I_off+(0.008*(1e6/si)))/(1e6/si)];
        for w = 2:length(k)
            if ((k(w) - k(w-1)) > double_counting) 
                eap_is = [eap_is, k(w)];
                large_yw = [large_yw, data(k(w) + I_off+(0.008*(1e6/si)), v_chan, sweep)];
                large_xw = [large_xw, (k(w) + I_off+(0.008*(1e6/si)))/(1e6/si)];
            end
        end
    end    
    eap_i{sweep} = eap_is; % get indices of spike onset (point where dV/dt exceeds 10 mV/ms)
    large_x{sweep} = large_xw;
    large_y{sweep} = large_yw;
    
    I = data(int16((I_on+I_off)/2),i_chan,sweep); % get current value for sweep
    I_hold = data(int16(I_on-100),i_chan,sweep);
    
    eap_table(sweep,1) = sweep; % record sweep #
    eap_table(sweep,2) = I; % record absolute current for sweep
    eap_table(sweep,3) = I - I_hold; % record current difference from hold for sweep
    eap_table(sweep,4) = size(eap_i{sweep},2); % record number of spikes for a sweep
end
clear v dv dv_dv k_v k_dv k eap_is I I_hold sweep w

eap_sweeps = eap_table(eap_table(:,4) > 0);
%% EAP properties

% Only run if there are EAPs
if isempty(eap_sweeps) == false
    eapPresence = "Yes";

    % Max number of EAPs & the corresponding sweep
    [maxNumEAPs, maxEAPSweep] = max(eap_table(:,4));

    % EAP threshold, or the number of evoked AIS APs required to elicit the
    % first EAP. Including APs elicited in ramp
    rampPath = append(extractBefore(filePath, "_eap.abf"), "_ramp.abf");
    [rampData,~,~]=abfload(rampPath);
    rampData = rampData(:,1); % I don't care about the current injection in this file
    % Detect APs in ramp
    db_rampData = gradient(rampData) ./ dt;
    k_vRamp = find(rampData > v_threshold);
    k_dvRamp = find(d_rampData > dv_dt_threshold);
    kRamp = intersect(k_vRamp, k_dvRamp);
    rampSpikeis = [kRamp(1)];
    for w = 2:length(kRamp)
        if (kRamp(w) - kRamp(w-1)) > double_counting
            rampSpikeis = [rampSpikeis, kRamp(w)];
        end
    end
    numRampSpikes = length(rampSpikeis);
    clear rampPath rampData d_rampData k_vRamp k_dvRamp kRamp rampSpikeis w
    eapThreshSweep = eap_sweeps(1);
    eapThresh = sum(spike_table(rheo_sweep:eapThreshSweep,4)) + numRampSpikes;
    % clear numRampSpikes

    % Make tables of EAPs, EAP ISIs, and barrages
    eap_props = zeros(0, 7); % Cols: sweep, index, time, amplitude, half-width, upstroke velocity, downstroke velocity
    isi_table = zeros(0, 3); % Cols: sweep, ISI (idx), ISI (ms)
    barrage_table = zeros(0,2); % Cols: sweep, barrage yes/no
    for i = 1:length(eap_sweeps) % for each sweep with EAPs
        sweep = eap_sweeps(i);
        for j = 1:eap_table(sweep,4) % for each EAP within the sweep
            idx = I_off+(0.008*(1e6/si)) + eap_i{sweep}(j);

            % Find amplitude by defining eap start time (idx) and eap end
            % time
            if j == eap_table(sweep,4)
                eapV = data(idx-0.002*(1e6/si):end,v_chan,sweep);
            else 
                otherside = 0.005*(1e6./si);
                eapV = data(idx-0.002*(1e6/si):(idx + otherside),v_chan,sweep);
            end
            clear otherside
            eapAmp = max(eapV) - data(idx-0.002*(1e6/si), v_chan, sweep);            

            % Calculate upstroke and downstroke velocity
            [~, eapPeaki] = max(eapV);
            eapVd = gradient(eapV)./(si*1e-3);
            eapUpVel = max(eapVd(1:eapPeaki));
            eapDownVel = min(eapVd(eapPeaki:end));
            clear eapVd 

            % Calculate EAP half-width
            halfUp = find(eapV >= (eapV(1) + 0.5*eapAmp));
            halfUp = halfUp(1);
            if (length(eapV) - eapPeaki) > 10
                halfDown = find(eapV(eapPeaki:end) <= eapV(1) + 0.5*eapAmp);
                halfDown = halfDown(1) + eapPeaki;
                eapHalfWidth = 1000 * (halfDown - halfUp) / (1e6/si); % units ms
            else
                eapHalfWidth = NaN;
            end

            % Populate properties table
            newRow = [sweep, idx, idx / (1e6/si), eapAmp, eapHalfWidth, eapUpVel, eapDownVel]; %, eapClass];
            eap_props = vertcat(eap_props, newRow);
            clear newRow eapAmp

            % Populate ISI table if there is more than one EAP in the sweep
            if j > 1
                newRowISI = [sweep, idx - (I_off + 100 + eap_i{sweep}(j-1)), (idx - (I_off + 100 + eap_i{sweep}(j-1)))/ (1e6/si)];
                isi_table = vertcat(isi_table, newRowISI);
            end
            clear newRowISI idx j
        end

        % Identify barrages (start within 50 ms, last 250 ms, 35 Hz) **
        % changed from 50 Hz
        meanEAPFreq = 1 / mean(isi_table(isi_table(:,1) == sweep, 3));
        if (eap_i{sweep}(1) < 0.05*1e6/si) && (eap_i{sweep}(end) > 0.25*1e6/si) && meanEAPFreq > 20
           barrage_table = vertcat(barrage_table, [sweep, "Barrage"]);
        else
           barrage_table = vertcat(barrage_table, [sweep, "None"]);
        end

    end
    clear i sweep

    % Calculate mean and SD ISI, if there are any
    if height(isi_table) > 0
        meanISI = mean(isi_table(:,3));
        sdISI = std(isi_table(:,3));
    else
        meanISI = "No sweeps with >1 EAP";
        sdISI = "No sweeps with >1 EAP";
    end

    % Detect barrage presence
    if ismember("Barrage", barrage_table(:,2)) == true
        barragePresence = "Yes";
    else
        barragePresence = "No";
    end
    clear barrage_table

else 
    eapPresence = "No";
    barragePresence = "No";
end