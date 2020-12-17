%load packet_1_QPSK.mat

%% Rx processing params
mod_type = 4;
%raw_rx_dec = data;
rx_data = raw_rx_dec;          % run OFDM tx code to get raw_rx_dec
LTS_CORR_THRESH = 0.8;         % Normalized threshold for LTS correlation
% Usage: Find all peaks whose magnitude is greater than 0.8 times
% the maximum magnitude after cross correlation (Packet Detection)

% Repeat the following code for each packet

%% Packet Detection
lts_xcorr = abs(xcorr(lts_t,sign(rx_data)));
lts_xcorr = fliplr(lts_xcorr(32:(length(rx_data)+32)));
max_lts_xcorr = max(lts_xcorr)*0.8;
lts_xcorr_peaks = find(lts_xcorr(1:end) > max_lts_xcorr);

payload_ind = max(lts_xcorr_peaks)+32;
lts_index = payload_ind-160;

% Output: Single packet extracted from rx_data
% with knowledge of preamble (LTS) indices and payload vector indices


%% CFO estimation and correction
% Use two copies of LTS for cross-correlation (Reference: Thesis)
FFT_OFFSET = 4;
rx_lts = rx_data(lts_index : lts_index+159);
rx_lts_1 = rx_data(lts_index+96 : lts_index+159);
rx_lts_2 = rx_data(lts_index+32 : lts_index+95);

N=64;
for n=1:N
    deltaf(n+1) = imag(rx_lts_1(n)/rx_lts_2(n));
    deltaf(n+1) = deltaf(n+1)/(2*pi*N);  
end 

rx_cfo_est_lts = mean(deltaf);     
rx_cfo_corr_t = exp(-1i*2*pi*rx_cfo_est_lts*[1:length(rx_data)]);
rx_dec_cfo_corr = rx_data .* rx_cfo_corr_t;
    

% Output: Packet with each value multiplied by CFO correction factor


%% CP Removal
% Refer to the process used to add CP at TX
% Converting vector back to matrix form will help

% Output: CP free payload matrix of size (N_SC * N_OFDM_SYMS)


%% FFT
% Refer to IFFT perfomed at TX
payload_rx = rx_dec_cfo_corr(payload_ind : payload_ind+N_OFDM_SYMS*(N_SC+CP_LEN)-1);
payload_matrix_form = reshape(payload_rx, (N_SC+CP_LEN), N_OFDM_SYMS);

% Remove the cyclic prefix, keeping FFT_OFFSET samples of CP (on average)
payload_mat_noCP = payload_matrix_form(CP_LEN - FFT_OFFSET +[1:N_SC], :);

% Take the FFT
syms_f_mat = fft(payload_mat_noCP, N_SC, 1);
% Output: Symbol matrix in frequency domain of same size


%% Channel estimation and correction
% Use the two copies of LTS and find channel estimate (Reference: Thesis)
% Convert channel estimate to matrix form and equlaize the above matrix
rx_lts = rx_dec_cfo_corr(lts_index : lts_index+159);
rx_lts_1 = raw_rx_dec(lts_index+96 : lts_index+159);
rx_lts_2 = raw_rx_dec(lts_index+32 : lts_index+95);

rx_lts1_f = fft(rx_lts_1);
rx_lts2_f = fft(rx_lts_2);
rx_lts_f = fft(rx_lts(-FFT_OFFSET + [97:160]));
rx_H_est = lts_f .* (rx_lts_f);

channel_matrix = repmat(transpose(rx_H_est),1,N_OFDM_SYMS);

syms_eq_mat = rdivide(syms_f_mat,channel_matrix); 

%syms_eq_mat = syms_f_mat ./ repmat(rx_H_est.', 1, N_OFDM_SYMS);
% Output : Symbol equalized matrix in frequency domain of same size


%% SFO estimation and correction using pilots
% SFO manifests as a frequency-dependent phase whose slope increases
% over time as the Tx and Rx sample streams drift apart from one
% another. To correct for this effect, we calculate this phase slope at
% each OFDM symbol using the pilot tones and use this slope to
% interpolate a phase correction for each data-bearing subcarrier.
pilots_f_mat = syms_eq_mat(SC_IND_PILOTS, :);
pilots_f_mat_comp = times(pilots_f_mat,pilots_mat);
    
% Calculate the phases of every Rx pilot tone
phase_angles = angle(pilots_f_mat_comp);
%pilot_phases = unwrap(angle(fftshift(pilots_f_mat_comp,1)), [], 1);

% Calculate slope of pilot tone phases vs frequency in each OFDM symbol
pilot_spacing_mat = repmat([14 50 14].', 1, N_OFDM_SYMS);
%pilot_spacing_
pilot_slope_mat = mean(diff(phase_angles) ./ pilot_spacing_mat);

% Calculate the SFO correction phases for each OFDM symbol
pilot_phase_sfo_corr = (1:64).' * pilot_slope_mat;
%pilot_phase_sfo_corr = repmat(transpose(pilot_slope_mat),1,N_OFDM_SYMS)
pilot_phase_corr = exp(-1i*(pilot_phase_sfo_corr));

% Apply the pilot phase correction per symbol
syms_eq_mat = syms_eq_mat .* pilot_phase_corr;
% Output: Symbol equalized matrix with pilot phase correction applied


%% Phase Error Correction using pilots
% Extract the pilots and calculate per-symbol phase error
pilots_f_mat = syms_eq_mat(SC_IND_PILOTS, :);
pilots_f_mat_comp = pilots_f_mat.*pilots_mat;
pilot_phase_err = angle(mean(pilots_f_mat_comp));
    
pilot_phase_err_corr = repmat(pilot_phase_err, N_SC, 1);
pilot_phase_corr = exp(-1i*(pilot_phase_err_corr));

% Apply the pilot phase correction per symbol
syms_eq_pc_mat = syms_eq_mat .* pilot_phase_corr;
payload_syms_mat = syms_eq_pc_mat(SC_IND_DATA, :);
rx_syms = reshape(payload_syms_mat, 1, N_DATA_SYMS);
% Output: Symbol equalized matrix with pilot phase correction applied
% Remove pilots and flatten the matrix to a vector rx_syms


%% Demodulation
figure(4);
scatter(real(rx_syms), imag(rx_syms),'filled');
title(' Signal Space of received bits');
xlabel('I'); ylabel('Q');
% FEC decoder
Demap_out = demapper(rx_syms,mod_type,1);

% viterbi decoder
rx_data_final = vitdec(Demap_out,trel,7,'trunc','hard');

% rx_data is the final output corresponding to tx_data, which can be used
% to calculate BER

[number,ber] = biterr(tx_data,rx_data_final);
ber
bit_errs = length(find(dec2bin(bitxor(tx_data, rx_data_final),8) == '1'));
bit_errs
