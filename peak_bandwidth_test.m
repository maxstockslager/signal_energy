function peak_bandwidth_test

clear all, close all

% plot settings
plot_peak_shapes = false;

% peak settings
peak_settings.transit_time = 10/1000; % seconds to transit the cantilever. transit time
peak_settings.mode_number = 1; % mode number. currently works for 1-3, could do higher if you know the eigenvalues
peak_settings.max_channel_length = 0.47; % cell travels this fraction of the cantilever length. 
    % set to ~0.95 for normal devices, or ~0.79 for the devices that stop around the node.
peak_settings.peak_height = 1; 
peak_settings.normalization = 'absolute'; % 'absolute' or 'antinodes';
    
% filter settings 
filter_settings.datarate = 5000; 
filter_settings.loop_order = 2; 
filter_settings.type = 'butterworth';

% set up sweep
bandwidth_sweep_range_Hz = 10:10:1000; 
signal_energy = zeros(size(bandwidth_sweep_range_Hz));
unfiltered_signal_energy = zeros(size(bandwidth_sweep_range_Hz));

% get the peak shape for this peak with these filter settings
for ii = 1 : numel(bandwidth_sweep_range_Hz)
    filter_settings.bandwidth = bandwidth_sweep_range_Hz(ii)*2*pi;
 
    peak_signals = get_peak_signals(peak_settings, filter_settings);
    if plot_peak_shapes, plot_signals(peak_signals), end
    signal_energy(ii) = calculate_signal_energy(peak_signals.time_full, ...
        peak_signals.filtSignal_full);
    unfiltered_signal_energy(ii) = calculate_signal_energy(peak_signals.time_full, ...
        peak_signals.idealSignal_full);

end

% get unfiltered signal energy to normalize
norm_signal_energy = signal_energy ./ unfiltered_signal_energy; 

% plot sweep results
figure
subplot(1, 2, 1)
plot(bandwidth_sweep_range_Hz, norm_signal_energy, 'k');
xlabel('Bandwidth (Hz)');
ylabel('Energy recovered');

subplot(1, 2, 2)
semilogy(bandwidth_sweep_range_Hz, norm_signal_energy, 'k');
xlabel('Bandwidth (Hz)');
ylabel('Energy recovered');

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function plot_signals(peak_signals)

figure
plot(peak_signals.time_vector*1000, peak_signals.df, 'k', 'LineWidth', 1.5)
hold on
xlabel('Time (ms)')
ylabel('Frequency shift (Hz)');
set(gca, 'YLim', [-1.05 0.05]*max(abs(get(gca, 'YLim'))));
set(gca, 'XLim', [0, max(peak_signals.time_vector)*1000]);

figure
plot([min(peak_signals.time_vector), max(peak_signals.time_vector)], [0 0], ...
    'Color', 0.65*[1 1 1]);
hold on
plot(peak_signals.time_vector, peak_signals.df, 'Color', 'k');
hold on
plot(peak_signals.time_vector, peak_signals.df_filt, 'r--');
% plot(peak_signals.timeVector_downsampled, peak_signals.df_downsampled, 'ro');
% legend('Frequency signal', '...with specified bandwidth', '...sampled at specified datarate');
xlabel('Time (s)');
ylabel('Relative frequency shift (Hz)');
% plot(get(gca, 'XLim'), [-1 -1], 'k--');
set(gca, 'YLim', [-1.05 0.05]);
% set(gcf, 'Position', [300 450 450 250])
set(gca, 'XLim', [min(peak_signals.time_vector), max(peak_signals.time_vector)]);
end

function energy = calculate_signal_energy(t, x)
energy = sum((abs(x)).^2);
end