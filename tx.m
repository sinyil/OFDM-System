function [txsignal, conf] = tx(txbits,conf,k)
% Digital Transmitter
%
%   [txsignal conf] = tx(txbits,conf,k) implements a complete transmitter
%   consisting of:
%       - modulator
%       - pulse shaping filter
%       - up converter
%   in digital domain.
%
%   txbits  : Information bits
%   conf    : Universal configuration structure
%   k       : Frame index
%

%% Mapping to BPSK or QPSK
if conf.modulation_order == 2
    GrayMap = 1/sqrt(2) * [(-1-1j), (-1+1j), (1-1j), (1+1j)];
    dualedBits = [txbits(1:2:end-1), txbits(2:2:end)]; % more convenient to map and demap.
    mapped = GrayMap(bi2de(dualedBits)+1).'; % to map
else
    error('Please use QPSK modulation.');
end

%% Training
training     = gentraining(conf.ntraining);
bpsktraining = -2*(training-0.5);

% Concatenate mapped signal with training
if conf.phaseTrack == 1 || conf.phaseTrack == 2
    trainedSignal = [bpsktraining; mapped];
elseif conf.phaseTrack == 3
    trainedSignal = double.empty();
    for i = 1:conf.nBlocks/2
         a = [bpsktraining; mapped((i-1)*conf.nsubchannels+1:i*conf.nsubchannels)];
         trainedSignal = [trainedSignal; a];
    end
elseif conf.phaseTrack == 4
    % For this type of tracking, we performed the parallellization here.
    % Training symbol is added later on, in parallelization block.
    untrrained_parallelForm = zeros(conf.nsubchannels,conf.nBlocks);
    pilot_seq = gentraining(conf.nPilots*conf.nBlocks);
    pilot_seq = -2*(pilot_seq-0.5);
    pilot_block = reshape(pilot_seq,[],conf.nBlocks);
    parallel_map = reshape(mapped, (conf.nsubchannels-conf.nPilots), []);
    % delta will be an even number.
    delta= 2*conf.nsubchannels/conf.nPilots;
    for i=1:conf.nPilots/2
        untrrained_parallelForm((i-1)*delta+1,:) = pilot_block(2*i-1,:);
        untrrained_parallelForm((i-1)*delta+2:i*delta-1,:) = parallel_map((i-1)*(delta-2)+1:i*(delta-2),:);
        untrrained_parallelForm(i*delta,:) = pilot_block(2*i,:);
    end
end

%% Parallelization
if conf.phaseTrack == 1 || conf.phaseTrack == 2 || conf.phaseTrack == 3 
   parallelForm = reshape(trainedSignal, conf.nsubchannels, []);
elseif conf.phaseTrack == 4
   parallelForm = [bpsktraining untrrained_parallelForm];
end

%% IFFT & Cyclic Prefix
[row, column] = size(parallelForm);
s = zeros(row*conf.os_factorOFDM, column);
cprefixed = zeros(conf.Ntotal*conf.os_factorOFDM, column);
for i = 1:column
    s(:,i) = osifft(parallelForm(:,i), conf.os_factorOFDM);
    cprefixed(:,i) = [s((row-conf.Ncp)*conf.os_factorOFDM+1:end,i); s(:,i)];
end

%% Serialization
serialForm = reshape(cprefixed, [], 1);

%% Preamble
preamble     = genpreamble(conf.npreamble);
bpskpreamble = -2*(preamble-0.5);

% Filtering preamble
upsampled_preamble = upsample(bpskpreamble, conf.os_factorPreamble);

rolloff_factor = 0.22;
filter = rrc(conf.os_factorPreamble, rolloff_factor, 20*conf.os_factorPreamble); %% remark : length of rrc
filtered_preamble = conv(upsampled_preamble, filter, 'same');

%% Energy Equalization between preamble and OFDM symbol
% To check if the number of the original OFDM symbols is smaller than 4
if column <= 4
symbols = serialForm(conf.Ntotal*conf.os_factorOFDM+1:end);
powerOFDM = sum(abs(symbols).^2)/length(symbols);
%timeOFDM = length(symbols)/conf.f_s;
else
% Take first three OFDM symbols coming after training symbol.
first4Symbols = serialForm(conf.Ntotal*conf.os_factorOFDM+1:5*conf.Ntotal*conf.os_factorOFDM);
powerOFDM = sum(abs(first4Symbols).^2)/length(first4Symbols);
%timeOFDM = 4*conf.nsubchannels*conf.os_factorOFDM/conf.f_s;
end

powerPreamble = sum(abs(filtered_preamble).^2)/length(filtered_preamble);
%timePreamble = length(filtered_preamble)/conf.f_s;

%A = sqrt(timePreamble*powerPreamble/timeOFDM/powerOFDM);
A = sqrt(powerPreamble/powerOFDM);

serialFormEnergyEqualized = serialForm * A;

%% Concatenate the preamble
preambledSignal = [filtered_preamble; serialFormEnergyEqualized];

%% Up Conversion
% to up convert and have real and imaginary parts of the signal in
% orthogonal spaces.
time = 1:1/conf.f_s:(length(preambledSignal)-1)/conf.f_s+1;
txsignal = real(preambledSignal .* exp(1i*2*pi*conf.f_c * time.'));
