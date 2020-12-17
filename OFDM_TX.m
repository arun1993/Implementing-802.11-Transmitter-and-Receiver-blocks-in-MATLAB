
clear;

%% Params:

CFO_FLAG = 1; % flag to enable CFO 
DETECTION_OFFSET = 100; % to add packet detection error

% Waveform params
N_OFDM_SYMS             = 500;         % Number of OFDM symbols
MOD_ORDER               =  4;          % Modulation order in power of 2 (1/2/4/6 = BSPK/QPSK/16-QAM/64-QAM)
TX_SCALE                = 1.0;         % Scale for Tx waveform ([0:1])

% OFDM params
SC_IND_PILOTS           = [8 22 44 58];                           % Pilot subcarrier indices
SC_IND_DATA             = [2:7 9:21 23:27 39:43 45:57 59:64];     % Data subcarrier indices
N_SC                    = 64;                                     % Number of subcarriers
CP_LEN                  = 16;                                     % Cyclic prefix length
N_DATA_SYMS             = N_OFDM_SYMS * length(SC_IND_DATA);      % Number of data symbols (one per data-bearing subcarrier per OFDM symbol)

channel_coding = .5; % coding rate 
trellis_end_length = 8; % bits for trellis to end 


%% Preamble
% Preamble is a concatenation of multiple copies of STS and LTS
% It is used for packet detection and CFO and channel estimation
% LTS is sufficient to be used for the above three blocks in a way similar to what is given in OFDM thesis.
% If you want to use STS in place of LTS, read the paper below:
% 'Robust Frequency and Timing Synchronization for OFDM' by Timothy M. Schmidl and Donald C. Cox

% STS
sts_f = zeros(1,64);
sts_f(1:27) = [0 0 0 0 -1-1i 0 0 0 -1-1i 0 0 0 1+1i 0 0 0 1+1i 0 0 0 1+1i 0 0 0 1+1i 0 0];
sts_f(39:64) = [0 0 1+1i 0 0 0 -1-1i 0 0 0 1+1i 0 0 0 -1-1i 0 0 0 -1-1i 0 0 0 1+1i 0 0 0];
sts_t = ifft(sqrt(13/6).*sts_f, 64);
sts_t = sts_t(1:16);

% LTS
lts_f = [0 1 -1 -1 1 1 -1 1 -1 1 -1 -1 -1 -1 -1 1 1 -1 -1 1 -1 1 -1 1 1 1 1 0 0 0 0 0 0 0 0 0 0 0 1 1 -1 -1 1 1 -1 1 -1 1 1 1 1 1 1 -1 -1 1 1 -1 1 -1 1 1 1 1];
lts_t = ifft(lts_f, 64);

preamble = [repmat(sts_t, 1, 30)  lts_t(33:64) lts_t lts_t];


%% Generate a payload of random integers
number_of_bits= (N_DATA_SYMS * MOD_ORDER - 2*trellis_end_length) * channel_coding;
tx_data = randi(2, 1, number_of_bits) - 1; 

% Forward Error Correction 
tx_data = double([tx_data zeros(1,trellis_end_length) ]);    % 8 bits padding
trel = poly2trellis(7, [171 133]);              % Define trellis
tx_code = convenc(tx_data,trel);            % convultional encoder 

% bits to signal space mapping these are you are x_k from the class
tx_syms = mapping(tx_code', MOD_ORDER, 1);

figure(1);
scatter(real(tx_syms), imag(tx_syms),'filled');
title(' Signal Space of transmitted bits');
xlabel('I'); ylabel('Q');

% Reshape the symbol vector to a matrix with one column per OFDM symbol,
tx_syms_mat = reshape(tx_syms, length(SC_IND_DATA), N_OFDM_SYMS);

% Define the pilot tone values as BPSK symbols
pilots = [1 1 -1 1].';

% Repeat the pilots across all OFDM symbols
pilots_mat = repmat(pilots, 1, N_OFDM_SYMS);


%% IFFT

% Construct the IFFT input matrix
ifft_in_mat = zeros(N_SC, N_OFDM_SYMS);

% Insert the data and pilot values; other subcarriers will remain at 0
ifft_in_mat(SC_IND_DATA, :)   = tx_syms_mat;
ifft_in_mat(SC_IND_PILOTS, :) = pilots_mat;

%Perform the IFFT --> frequency to time translation
tx_payload_mat = ifft(ifft_in_mat, N_SC, 1);

% Insert the cyclic prefix
if(CP_LEN > 0)
    tx_cp = tx_payload_mat((end-CP_LEN+1 : end), :);
    tx_payload_mat = [tx_cp; tx_payload_mat];
end

% Reshape to a vector
tx_payload_vec = reshape(tx_payload_mat, 1, numel(tx_payload_mat));

% Construct the full time-domain OFDM waveform
tx_vec = [preamble tx_payload_vec];


%% Interpolate, 
% Interpolation filter basically implements the DAC before transmission
% On the receiver's end decimation is performed to implement the ADC

% Define a half-band 2x interpolation filter response
interp_filt2 = zeros(1,43);
interp_filt2([1 3 5 7 9 11 13 15 17 19 21]) = [12 -32 72 -140 252 -422 682 -1086 1778 -3284 10364];
interp_filt2([23 25 27 29 31 33 35 37 39 41 43]) = interp_filt2(fliplr([1 3 5 7 9 11 13 15 17 19 21]));
interp_filt2(22) = 16384;
interp_filt2 = interp_filt2./max(abs(interp_filt2));

% Pad with zeros for transmission to deal with delay through the interpolation filter
tx_vec_padded = [tx_vec, zeros(1, ceil(length(interp_filt2)/2))];
tx_vec_2x = zeros(1, 2*numel(tx_vec_padded));
tx_vec_2x(1:2:end) = tx_vec_padded;
tx_vec_air = filter(interp_filt2, 1, tx_vec_2x);

figure(2);
plot(abs(tx_vec_2x));
hold on;
plot(abs(tx_vec_air(22:end)));
xlim([20,50]);
title('Interpolation visualized');
xlabel('time'); ylabel('amplitude');
legend('y = pre filtering','y = post filtering')

% Scale the Tx vector to +/- 1, becasue ADC and DAC take samples input from
% 1 to -1
tx_vec_air = TX_SCALE .* tx_vec_air ./ max(abs(tx_vec_air));

figure(3);
plot(db(abs(fftshift(fft(tx_vec_air)))));
xlim([20000,60000]); ylim([0,65]);
% in this plot, why do see four peaks?


%% This part of code is for simulating the wireless channel.
% You can later use the receiver raw data file instead to test your code.

% Perfect (ie. Rx=Tx):
% rx_vec_air = tx_vec_air;

% to enable CFO make CFO_FLAG=1
if(CFO_FLAG)
    tx_vec_air = tx_vec_air .* exp(-1i*2*pi*1e-4*[0:length(tx_vec_air)-1]);
end

% AWGN:

TRIGGER_OFFSET_TOL_NS   = 3000;        % Trigger time offset toleration between Tx and Rx that can be accomodated
SAMP_FREQ               = 40e6;        % Sampling frequency

rx_vec_air = [tx_vec_air, zeros(1,ceil((TRIGGER_OFFSET_TOL_NS*1e-9) / (1/SAMP_FREQ)))];
rx_vec_air = rx_vec_air + 0*complex(randn(1,length(rx_vec_air)), randn(1,length(rx_vec_air)));


% Decimate
raw_rx_dec = filter(interp_filt2, 1, rx_vec_air);
raw_rx_dec = [zeros(1,DETECTION_OFFSET) raw_rx_dec(1:2:end)];


