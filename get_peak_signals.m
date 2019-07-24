function peak_signals = get_peak_signals(peak_settings, filter_settings)

transit_time = peak_settings.transit_time;
peak_height = peak_settings.peak_height; 
max_channel_length = peak_settings.max_channel_length; 

n = peak_settings.mode_number; 

datarate = filter_settings.datarate; 

% get the mode shape and peak shape 
[u, ~] = get_mode_shape(n, max_channel_length);
df = get_peak_shape(u);
df = df * peak_height; 
dt = peak_settings.transit_time/length(df); % still correct after change
timeToPadZeros = max([4*1/filter_settings.bandwidth, peak_settings.transit_time]); % sec 
zerosToPad = ceil(timeToPadZeros/dt); 

df_padded = [zeros(1, zerosToPad), df, zeros(1, zerosToPad)];
[time_vector, dt] = get_time_vector(transit_time, numel(df), numel(df_padded));
df_filt = filter_signal(time_vector, df_padded, filter_settings);
[t_downsampled, signal_downsampled] = random_downsample(time_vector, df_filt, datarate);

peak_signals.time_full = time_vector;
peak_signals.idealSignal_full = df_padded; 
peak_signals.filtSignal_full = df_filt; 
peak_signals.time_downsampled = t_downsampled;
peak_signals.filt_signal_downsampled = signal_downsampled; 
peak_signals.dt_full = dt;
peak_signals.dt_ds = 1/datarate; 

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [u, x] = get_mode_shape(n, max_channel_length)
% n is the mode number (currently works for 1-3)
% max_channel_length is the fraction of the cantilever length that the cell travels before
% turning around (0 < L < 1)
% evaluated at the grid in x. 

lambda = [1.8751, 4.6941, 7.855]; % eigenvalues 
L = 1; % cantilever length
An = L/1000; % amplitude. doesn't matter
dx = 0.001; % step size to evaluate mode shape
x = 0:dx:max_channel_length*L; % vector to evaluate mode shape 
x = [x, x(end:-1:1)]; % add points going the opposite direction
u = An/2 * ( (cosh(lambda(n)*x/L) - cos(lambda(n)*x/L)) - ((cosh(lambda(n))+cos(lambda(n)))/(sinh(lambda(n))+sin(lambda(n)))) * ...
    (sinh(lambda(n)*x/L) - sin(lambda(n)*x/L))); % mode shape
u = u/max(u); % normalize to max amplitude of 1 
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function df = get_peak_shape(u)
% get the peak shape resulting from a cell traveling through the cantilever
% with the mode shape given by u. normalized to 1. 

% what's the resonance frequency shift resulting from the added mass?
dm = 0.01; % added mass. magnitude doesn't actually matter since we just want the peak shape. 
m = 1; % cantilever mass 
df = -1 + (1 + u.^2 * dm/m).^-0.5;

% optional: normalize so that antinodes have height 1
% df = df/abs(min((df(1 : round(length(df)/3)))));
df = df/abs(min(df));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [time_vector, dt] = get_time_vector(transit_time, signal_length, padded_length)
% creates time vector to simulate the particle going through the cantilever
% with the specified transit time.
% signal_length is the length of the *peak*. 
% padded_length is the total length of the signal vector, i.e., padded with
% zeros on either side. 
idx_padded = 0 : (padded_length-1);
dt = transit_time/signal_length;
time_vector = idx_padded * dt;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [t_downsampled, signal_downsampled] = random_downsample(t, signal, datarate)
% downsamples the signal (evaluated at time vector t) at the given datarate. 
% random phase between signal and sampling. 
timeStepSize = (t(end) - t(1)) / numel(t);
sampleTime = 1/datarate;
downsampleRatio = ceil(sampleTime / timeStepSize);
tempRand = rand();
signal_downsampled = signal(ceil(downsampleRatio*tempRand) : downsampleRatio : end);
t_downsampled = t(ceil(downsampleRatio*tempRand) : downsampleRatio : end);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function signal_filtered = filter_signal(t, signal, filter_settings)
% input bandwidth in rad/s (NOT Hz)
s = tf('s');

switch filter_settings.type
    case 'butterworth'     
        [z, p, k] = butter(filter_settings.loop_order, filter_settings.bandwidth, 's');
        [num, den] = zp2tf(z, p, k);
        SMR_TF = tf(num,den);
    case 'chebyshev'
       [z,p,k] = cheb1ap(filter_settings.loop_order, filter_settings.ripple);       % Lowpass filter prototype
       z = z * filter_settings.bandwidth;
       p = p * filter_settings.bandwidth; 
       [num,den] = zp2tf(z,p,k);     % Convert to transfer function form 
       SMR_TF = tf(num, den);
end

signal_filtered = lsim(SMR_TF, signal, t);
end