
%% Rx processing params
%load packet_1_QPSK.mat
%raw_rx_dec = data;
MOD_ORDER = 2;
rx_data = raw_rx_dec;          % run OFDM tx code to get raw_rx_dec
LTS_CORR_THRESH = 0.8;         % Normalized threshold for LTS correlation
% Usage: Find all peaks whose magnitude is greater than 0.8 times
% the maximum magnitude after cross correlation (Packet Detection)

% Repeat the following code for each packet

%% Packet Detection
% Cross correlation of received signal with LTS
lts_xcorr = abs(xcorr(lts_t,sign(rx_data)));
%lts_xcorr = fliplr(lts_xcorr(32:(length(rx_data)+32)));
max_lts_xcorr = max(lts_xcorr)*0.8;
lts_xcorr_peaks = find(lts_xcorr(1:end) > max_lts_xcorr);

payload_index = max(lts_xcorr_peaks)+64;
lts_index = payload_index-160;
rx_data(payload_index:end)

% Output: Single packet extracted from rx_data
% with knowledge of preamble (LTS) indices and payload vector indices


%% CFO estimation and correction
% Use two copies of LTS for cross-correlation (Reference: Thesis)
rx_lts = rx_data(lts_index : lts_index+159);
rx_lts_1 = rx_data(lts_index+96 : lts_index+159);
rx_lts_2 = rx_data(lts_index+32 : lts_index+95);

N=64;
for n=1:N
    deltaf(n+1) = imag(rx_lts_1(n)/rx_lts_2(n));
    deltaf(n+1) = deltaf(n+1)/(2*pi*N);  
end 

rx_cfo_estimate= mean(deltaf);     
rx_cfo_correction = exp(-1i*2*pi*rx_cfo_estimate*[1:length(rx_data)]);
rx_data_cfo_corrected = times(rx_data, rx_cfo_correction);
rx_data_cfo_corrected(payload_index:end)
% Output: Packet with each value multiplied by CFO correction factor


%% CP Removal
% Refer to the process used to add CP at TX
% Converting vector back to matrix form will help

FFT_OFFSET = 4;
payload_received = rx_data_cfo_corrected(payload_index : payload_index + N_OFDM_SYMS * (N_SC + CP_LEN) - 1);
payload_matrix_form = reshape(payload_received, (N_SC+CP_LEN), N_OFDM_SYMS);
payload_remove_CP = payload_matrix_form(CP_LEN - FFT_OFFSET +[1:N_SC], :);
payload_remove_CP
% Output: CP free payload matrix of size (N_SC * N_OFDM_SYMS)


%% FFT
% Refer to IFFT perfomed at TX
symbol_matrix = fft(payload_remove_CP, N_SC, 1);
symbol_matrix
% Output: Symbol matrix in frequency domain of same size


%% Channel estimation and correction
% Use the two copies of LTS and find channel estimate (Reference: Thesis)
% Convert channel estrx_lts = rx_dec_cfo_corr(lts_index : lts_index+159);
rx_lts_1 = rx_data(lts_index+96 : lts_index+159);
rx_lts_2 = rx_data(lts_index+32 : lts_index+95);

rx_lts1_f = fft(rx_lts_1);
rx_lts2_f = fft(rx_lts_2);
rx_lts_f = fft(rx_lts(-FFT_OFFSET + [97:160]));
channel_estimate = times(lts_f , rx_lts_f);
channel_matrix = repmat(transpose(channel_estimate),1,N_OFDM_SYMS);
symbol_equalized_matrix = rdivide(symbol_matrix,channel_matrix); 

% Output : Symbol equalized matrix in frequency domain of same size


%% SFO estimation and correction using pilots
% SFO manifests as a frequency-dependent phase whose slope increases
% over time as the Tx and Rx sample streams drift apart from one
% another. To correct for this effect, we calculate this phase slope at
% each OFDM symbol using the pilot tones and use this slope to
% interpolate a phase correction for each data-bearing subcarrier.
pilots_f = symbol_equalized_matrix(SC_IND_PILOTS, :);
pilots_f_matrix = times(pilots_f,pilots_mat);
phase_angle_drift = angle(pilots_f_matrix);

% pilot_slope_matrix = diff(phase_angle_drift) ./ diff(pilots_f_matrix);

pilot_spacing_matrix = repmat(transpose([14 50 14]), 1, N_OFDM_SYMS);
pilot_slope_matrix = mean(diff(phase_angle_drift) ./ pilot_spacing_matrix);

pilot_symbols_representation = (1:64);
pilot_phase_sfo_correction = times(transpose(pilot_symbols_representation), pilot_slope_matrix);
pilot_phase_correction_matrix = exp(-1i*(pilot_phase_sfo_correction));
symbol_equalized_matrix = times(symbol_equalized_matrix, pilot_phase_correction_matrix);

symbol_equalized_matrix
% Output: Symbol equalized matrix with pilot phase correction applied


%% Phase Error Correction using pilots
% Extract the pilots and calculate per-symbol phase error
pilots_f = symbol_equalized_matrix(SC_IND_PILOTS, :);
pilots_f_matrix = times(pilots_f,pilots_mat);

pilot_phase_error = angle(mean(pilots_f_matrix));
    
pilot_phase_error_correction = repmat(pilot_phase_error, N_SC, 1);
pilot_phase_correction_matrix = exp(-1i*(pilot_phase_error_correction));

symbols_phasecorrection_matrix = times(symbol_equalized_matrix, pilot_phase_correction_matrix);
payload_symbols_matrix = symbols_phasecorrection_matrix(SC_IND_DATA, :);
rx_syms = reshape(payload_symbols_matrix, 1, N_DATA_SYMS);
rx_syms


% Output: Symbol equalized matrix with pilot phase correction applied
% Remove pilots and flatten the matrix to a vector rx_syms


%% Demodulation

figure(4);
scatter(real(rx_syms), imag(rx_syms),'filled');
title(' Signal Space of received bits');
xlabel('I'); ylabel('Q');

% FEC decoder
Demap_out = demapper(rx_syms,MOD_ORDER,1);

% viterbi decoder
rx_data_final = vitdec(Demap_out,trel,7,'trunc','hard');

% rx_data is the final output corresponding to tx_data, which can be used
% to calculate BER

[number,ber] = biterr(tx_data,rx_data_final);
ber

