function [rxbits, conf] = rx(rxsignal,conf,k)
% Digital Receiver
%
%   [txsignal conf] = tx(txbits,conf,k) implements a complete causal
%   receiver in digital domain.
%
%   rxsignal    : received signal
%   conf        : configuration structure
%   k           : frame index
%
%   Outputs
%
%   rxbits      : received bits
%   conf        : configuration structure
%

rxbits = zeros(conf.nbits,1);

%% Down Conversion
time = 1:1/conf.f_s:(length(rxsignal)-1)/conf.f_s+1;
downConverted = rxsignal .* exp(1i*-2*pi*(conf.f_c+conf.offset) * time.');

%% Low Pass Filtering
fCorner = ceil((conf.nsubchannels+1)/2)*conf.f_spacing;
lowPassed = 2 * ofdmlowpass(downConverted, conf, fCorner);

%% Matched Filter

rolloff_factor = 0.22;
filter = rrc(conf.os_factorPreamble, rolloff_factor, 20*conf.os_factorPreamble);
filtered_rx_signal = conv(lowPassed, filter, 'same');

%% Frame Synchronization
[data_idx] = frame_sync(filtered_rx_signal, conf.os_factorPreamble, conf.npreamble); % Index of the first data symbol

%% Parallelization
data = lowPassed(data_idx:data_idx+conf.os_factorOFDM*conf.Ntotal*(conf.nBlocks+1)-1);
parallelForm = reshape(data, conf.Ntotal*conf.os_factorOFDM, []);

%% Removing Cyclic Prefix & FFT
[~, column] = size(parallelForm);
downFreqData = zeros(conf.nsubchannels, column);
for i = 1:column
    CPRemoved = parallelForm(conf.Ncp*conf.os_factorOFDM+1:end,i);
    downFreqData(:,i) = osfft(CPRemoved, conf.os_factorOFDM);
end

%% Training
training     = gentraining(conf.ntraining);
bpsktraining = -2*(training-0.5);

%% Channel Estimation and Correction by using one of conf.phaseTrack methods.
[excluded, HHat, hhat, delaySpread] = phaseTrack_TrainingExclude(downFreqData, bpsktraining, conf);

%% Serialization
serialForm = reshape(excluded, [], 1);

%% Demapper
if conf.modulation_order == 2
    GrayMap = 1/sqrt(2) * [(-1-1j), (-1+1j), (1-1j), (1+1j)];
    [~,ind] = min((ones(conf.nbits/2,4)*diag(GrayMap) - serialForm(:,[1 1 1 1])),[],2);
    demapped = de2bi(ind-1);
    rxbits(1:2:end-1) = demapped(:,1);
    rxbits(2:2:end) = demapped(:,2);
end
end