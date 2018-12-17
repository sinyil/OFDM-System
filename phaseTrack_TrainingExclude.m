function [excluded, Hhat, hhat, delaySpread] = phaseTrack_TrainingExclude(downFreqData, bpsktraining, conf)

    [~, column] = size(downFreqData);
    
    if conf.phaseTrack == 1
        % Channel Estimation
        bpskTrainingDiag = diag(bpsktraining);
        Z = downFreqData(:,1);
        Hhat = inv(bpskTrainingDiag) * Z; % Channel Transfer Function
        hhat = ifft(Hhat); % Channel Impulse Response estimate duw to trainning
        delaySpread = find(diff(abs(hhat))>0,1,'first') * conf.os_factorOFDM / conf.f_s;

        % Plot the channel impulse response estimate
        figure; plot(conf.frequencyVector,abs(Hhat)); title('$$|\hat{H}|$$','Interpreter','Latex'); xlabel('Frequency (Hz)');
        figure; plot(conf.frequencyVector,angle(Hhat)); title('$$Angle of \hat{H}$$','Interpreter','Latex'); xlabel('Frequency (Hz)');
        figure; plot(abs(hhat)); title('$$|\hat{h}|$$','Interpreter','Latex');
        
        % Channel correction for training symbol
        correctedSignal = downFreqData ./ Hhat;
        
        % To exclude the training symbol located at the beginning.
        excluded = correctedSignal(:,2:end);
        
    elseif conf.phaseTrack == 2
        %% First Channel Estimation 
        bpskTrainingDiag = diag(bpsktraining);
        Z = downFreqData(:,1);
        Hhat = inv(bpskTrainingDiag) * Z; % Channel Transfer Function
        hhat = ifft(Hhat); % Channel Impulse Response estimate duw to trainning
        delaySpread = find(abs(hhat)<0.1,1,'first') * conf.os_factorOFDM / conf.f_s;

        % Plot the channel impulse response estimate
        figure; plot(conf.frequencyVector,abs(Hhat)); title('$$|\hat{H}|$$','Interpreter','Latex'); xlabel('Frequency (Hz)');
        figure; plot(conf.frequencyVector,angle(Hhat)); title('$$Angle of \hat{H}$$','Interpreter','Latex'); xlabel('Frequency (Hz)');
        figure; plot(abs(hhat)); title('$$|\hat{h}|$$','Interpreter','Latex');

        %% Channel Amplitude Correction
        ampCorrected = zeros(conf.nsubchannels,column);
        correctedSignal = zeros(conf.nsubchannels,column);

        % Channel correction for training symbol
        correctedSignal(:,1) = downFreqData(:,1) ./ Hhat;
        ampCorrected(:,1) = downFreqData(:,1) ./ abs(Hhat);

        for i = 2:column
            ampCorrected(:,i) = downFreqData(:,i) ./ abs(Hhat);
        end

        %% Channel Phase Correction and Phase Tracking with Viterbi-Viterbi
        theta_hat = zeros(conf.nsubchannels,column);
        theta_hat(:,1) = angle(Hhat);

        for i=2:column
          % Phase estimation    
          % Apply viterbi-viterbi algorithm
            deltaTheta = 1/4*angle(-ampCorrected(:,i).^4) + pi/2*(-1:4);

          % Unroll phase
            [~, ind] = min(abs(deltaTheta - theta_hat(:,i-1)),[],2);
            theta = deltaTheta(ind);

          % Lowpass filter phase
            theta_hat(:,i) = mod(0.01*theta + 0.99*theta_hat(:,i-1), 2*pi);

          % Phase correction
            correctedSignal(:,i) = ampCorrected(:,i) .* exp(-1j * theta_hat(:,i));   % ...and rotate the current symbol accordingly
        end
        
        figure;
        plot(theta_hat); title('$$\hat{\theta}$$','Interpreter','Latex'); xlabel('Frequency (Hz)');

        %% Exclude Training Symbol and Serialization
        excluded = correctedSignal(:,2:end);
    elseif conf.phaseTrack == 3
        % Channel Estimation
        bpskTrainingDiag = diag(bpsktraining);
        Hhat = zeros(conf.nsubchannels,conf.nBlocks/2);
        hhat = zeros(conf.nsubchannels,conf.nBlocks/2);
        delaySpread = zeros(conf.nBlocks/2);
        for i = 1:2:column-1
            Z = downFreqData(:,i);
            Hhat(:,(i+1)/2) = inv(bpskTrainingDiag) * Z; % Channel Transfer Function
            hhat(:,(i+1)/2) = ifft(Hhat(:,(i+1)/2)); % Channel Impulse Response estimate due to training
            delaySpread((i+1)/2) = find(diff(abs(hhat))>0,1,'first') * conf.os_factorOFDM / conf.f_s;
        end
        
        % Channel correction for training symbol
        ofdmSignals = downFreqData(:,2:2:end);
        correctedSignal = ofdmSignals ./ Hhat;
        
        % To exclude the training symbol located at the beginning.
        excluded = correctedSignal;
        
        figure;
        for i = 1:conf.nBlocks/2
            % Plot the amplitude of the channel transfer function
            plot(conf.frequencyVector,abs(Hhat));
            title('$$|\hat{H}|$$','Interpreter','Latex'); xlabel('Frequency (Hz)');
            hold on;
        end
        
        figure;
        for i = 1:conf.nBlocks/2
            % Plot the angle of the channel transfer function
            plot(conf.frequencyVector,angle(Hhat));
            title('$$Angle of \hat{H}$$','Interpreter','Latex'); xlabel('Frequency (Hz)');
            hold on;
        end
        
        figure;
        for i = 1:conf.nBlocks/2
            % Plot the channel impulse response
            plot(abs(hhat));
            title('$$|\hat{H}|$$','Interpreter','Latex');
            hold on;
        end
    elseif conf.phaseTrack == 4
        %% First Channel Estimation 
        bpskTrainingDiag = diag(bpsktraining);
        Z = downFreqData(:,1);
        Hhat = inv(bpskTrainingDiag) * Z; % Channel Transfer Function
        hhat = ifft(Hhat); % Channel Impulse Response estimate duw to trainning
        delaySpread = find(abs(hhat)<0.1,1,'first') * conf.os_factorOFDM / conf.f_s;
        
        % Plot the channel impulse response estimate
        figure; plot(abs(Hhat)); title('Hhat Abs');
        figure; plot(angle(Hhat)); title('Hhat Angle');
        figure; plot(abs(hhat)); title('hhat Abs');

        %% Channel Amplitude Correction
        ampCorrected = zeros(conf.nsubchannels,column);
        correctedSignal = zeros(conf.nsubchannels-conf.nPilots,column-1);
        for i = 1:column
            ampCorrected(:,i) = downFreqData(:,i)./Hhat;
        end
        
        %% Erasing Pilot Symbols and Correction
        org_dataFreq = zeros(conf.nsubchannels-conf.nPilots,column-1);
        % Also taken from pseudo random code.
        pilot = gentraining(conf.nPilots*(column-1));
        pilot = -2*(pilot-0.5);
        pilot_ref = reshape(pilot,[],column-1);
        pilot_data = zeros(conf.nPilots,column-1);
        delta= 2*conf.nsubchannels/conf.nPilots;
        
        % Reconstruction of received pilots and original data
        for i=1:conf.nPilots/2
            pilot_data(2*i-1,:) = ampCorrected((i-1)*delta+1,2:column);
            org_dataFreq((i-1)*(delta-2)+1:i*(delta-2),:) = ampCorrected((i-1)*delta+2:i*delta-1,2:column);
            pilot_data(2*i,:) = ampCorrected(i*delta,2:column);
            pilot_H1 = pilot_data(2*i-1,:)./pilot_ref(2*i-1,:);
            pilot_H2 = pilot_data(2*i,:)./pilot_ref(2*i,:) ;
            
            % To interpolate linearly and estimate the phase at the subcarriers that
            % are located between pilot subcarriers.
            for kk = 1:delta-2
                delta_theta = ((delta-kk).*angle(pilot_H1)-(kk.*angle(pilot_H2)))/(delta);
                correctedSignal(kk + (i-1)*(delta-2),:) = org_dataFreq(kk+(i-1)*(delta-2),:).*exp(-1j *delta_theta);
            end
        end
        % Training and pilot symbols are already excluded.
        excluded = correctedSignal;  
    end

end
