% % % % %
% Wireless Receivers: algorithms and architectures
% Audio Transmission Framework 
%
%
%   3 operating modes:
%   - 'matlab' : generic MATLAB audio routines (unreliable under Linux)
%   - 'native' : OS native audio system
%       - ALSA audio tools, most Linux distrubtions
%       - builtin WAV tools on Windows 
%   - 'bypass' : no audio transmission, takes txsignal as received signal

% Configuration Values
conf.audiosystem = 'matlab'; % Values: 'matlab','native','bypass'
conf.inputType = 'random'; % 'random', 'image': sends the image in 8 frames. (8 bits per pixel)

conf.f_s     = 48000;   % sampling rate
conf.f_sym   = 100;     % symbol rate
conf.nframes = 1;       % number of frames to transmit
conf.tOFDMBlock = 30; % Twice the number of OFDM blocks. (QPSK)
conf.modulation_order = 2; % QPSK:2
conf.f_c     = 8000; % carrier frequency
conf.nsubchannels = 512; % number of subcarriers. For phase tracking 4, # of subcarriers should be even.
conf.f_spacing   = 5; % frequency spacing between subcarriers.
conf.phaseTrack = 1; % Only Training: 1, Phase Tracking: 2, Training per OFDM Symbol = 3, Pilot Tracking = 4
conf.nPilots = 16 ; % For only Pilot Tracking. Compatible with even numbers only.

if strcmp(conf.inputType,'random')
    if conf.phaseTrack == 4
        conf.nbits   = (conf.nsubchannels-conf.nPilots)*conf.tOFDMBlock;    % number of bits per frame. Chosen to be multiple of # of subcarriers.
    else
        conf.nbits   = conf.nsubchannels*conf.tOFDMBlock;    % number of bits per frame. Chosen to be multiple of # of subcarriers.
    end 
elseif strcmp(conf.inputType,'image')
    load('lena64.mat');
    [r, c] = size(lena_64x64);
    binstring = dec2bin(lena_64x64, 8);
    txbitsAll = reshape(str2num(binstring(:)),r*c,8);
    [nb, nf] = size(txbitsAll);
    conf.nbits = nb;
    conf.nframes = nf;
end

conf.frequencyVector = conf.f_spacing:conf.f_spacing:conf.f_spacing*conf.nsubchannels;
conf.os_factorPreamble = conf.f_s / conf.f_sym;
conf.os_factorOFDM = conf.f_s / (conf.f_spacing * conf.nsubchannels);
conf.Ncp = 20; % minimum value of Ncp is 6. It is limited by channel impulse response.
conf.Ntotal = conf.nsubchannels + conf.Ncp;

if conf.phaseTrack == 1 || conf.phaseTrack == 2
    conf.nBlocks = (conf.nbits/conf.nsubchannels) / conf.modulation_order;
elseif conf.phaseTrack == 3
    conf.nBlocks = 2 * (conf.nbits/conf.nsubchannels) / conf.modulation_order;
elseif conf.phaseTrack == 4
    if strcmp(conf.inputType,'random')
        conf.nBlocks = conf.nbits/(conf.nsubchannels-conf.nPilots) / conf.modulation_order;
    elseif strcmp(conf.inputType,'image')
        conf.nBlocks = length(txbitsAll(1:end-(conf.nPilots*length(txbitsAll)/conf.nsubchannels),:)) / (conf.nsubchannels-conf.nPilots) / conf.modulation_order;
        conf.nbits = length(txbitsAll(1:end-(conf.nPilots*length(txbitsAll)/conf.nsubchannels),:));
    end
else
    error('No phase tracking method found. Enter an appropriate method.');
end

conf.ntraining = conf.nsubchannels;
conf.npreamble  = 100;
conf.bitsps     = 16;   % bits per audio sample
conf.offset     = 0;

% Init Section
% all calculations that you only have to do once
conf.os_factor = conf.f_s/conf.f_sym;
if mod(conf.os_factor,1) ~= 0
   disp('WARNING: Sampling rate must be a multiple of the symbol rate'); 
end
conf.nsyms      = ceil(conf.nbits/conf.modulation_order);

% Initialize result structure with zero
res.biterrors   = zeros(conf.nframes,1);
res.rxnbits     = zeros(conf.nframes,1);

% TODO: To speed up your simulation pregenerate data you can reuse
% beforehand.

% Results

for k=1:conf.nframes
    
    if strcmp(conf.inputType,'random')
        % Generate random data
        txbits = randi([0 1],conf.nbits,1);
    elseif strcmp(conf.inputType,'image')
        % Take the first column of the image data as frame.
        if conf.phaseTrack == 4
            % To get the compatible number of bits, exclude some of the last pixels.
            txbits = txbitsAll(1:end-(conf.nPilots*length(txbitsAll)/conf.nsubchannels),k);
        else
            txbits = txbitsAll(:,k);
        end
    end
    
    % TODO: Implement tx() Transmit Function
    [txsignal, conf] = tx(txbits,conf,k);
    
    % Plot transmitted signal for debugging
    figure;
    plot(txsignal);
    title('Transmitted Signal')
    
    % % % % % % % % % % % %
    % Begin
    % Audio Transmission
    %
    
    % normalize values
    peakvalue       = max(abs(txsignal));
    normtxsignal    = txsignal / (peakvalue + 0.3);
    
    % create vector for transmission
    rawtxsignal = [ zeros(conf.f_s,1) ; normtxsignal ;  zeros(conf.f_s,1) ]; % add padding before and after the signal
    rawtxsignal = [  rawtxsignal  zeros(size(rawtxsignal)) ]; % add second channel: no signal
    txdur       = length(rawtxsignal)/conf.f_s; % calculate length of transmitted signal
    
    %wavwrite(rawtxsignal,conf.f_s,16,'out.wav')   
    audiowrite('out.wav',rawtxsignal,conf.f_s)  
        
    % MATLAB audio mode
    if strcmp(conf.audiosystem,'matlab')
        disp('MATLAB generic');
        playobj = audioplayer(rawtxsignal,conf.f_s,conf.bitsps);
        recobj  = audiorecorder(conf.f_s,conf.bitsps,1);
        record(recobj);
        disp('Recording in Progress');
        playblocking(playobj)
        pause(0.5);
        stop(recobj);
        disp('Recording complete')
        rawrxsignal  = getaudiodata(recobj,'int16');
        rxsignal     = double(rawrxsignal(1:end))/double(intmax('int16')) ;
        
    elseif strcmp(conf.audiosystem,'bypass')
        rawrxsignal = rawtxsignal(:,1);
        rxsignal    = rawrxsignal;
    end
    
    % Plot received signal for debugging
    figure;
    plot(rxsignal);
    title('Received Signal')
    
    %
    % End
    % Audio Transmission   
    % % % % % % % % % % % %
    
    % TODO: Implement rx() Receive Function
    [rxbits, conf]       = rx(rxsignal,conf,k);
    
    if strcmp(conf.inputType,'image')
        eightBits(:,k) = rxbits;
    end
    
    res.rxnbits(k)      = length(rxbits);  
    res.biterrors(k)    = sum(rxbits ~= txbits);
    
end

per = sum(res.biterrors > 0)/conf.nframes;
ber = sum(res.biterrors)/sum(res.rxnbits);

if strcmp(conf.inputType,'image')
    % Create and display the reimagined.
    decimals = bi2de(eightBits, 'left-msb');
    if conf.phaseTrack == 4
        % Deleted Pixels will be recovered as white pixel.
        deletedPixels = 255*ones(conf.nPilots*length(txbitsAll)/conf.nsubchannels,1);
        reimagined = reshape([decimals; deletedPixels],r,c);
    else
        reimagined = reshape(decimals,r,c);
    end
    figure; imshow(uint8(reimagined)); title('Reimagined');
    pixelErrorRate = sum(sum(uint8(reimagined) ~= lena_64x64))/(r*c);
end
