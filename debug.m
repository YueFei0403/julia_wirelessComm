clc;
clear all;
close all;
warning('off')

%% === Configure Base-station transmitters  ===
carrier = nrCarrierConfig;
carrier.NSizeGrid = 1; % Bandwidth in RB
carrier.SubcarrierSpacing = 15;
carrier.CyclicPrefix = 'Normal';
ofdmInfo = nrOFDMInfo(carrier);

%% === Configure Transmit and Receive Antenna Arrays ===
simParams.fc = 6e9; % freqRange = 'FR1';
simParams.c = physconst('LightSpeed');

simParams.lambda = simParams.c/simParams.fc;
simParams.NUser = 1;
simParams.NumTx = 1; 
simParams.NumRx = 8;
simParams.NumPaths = 3;
simParams.NumDU = 4; % intensity of DUs -- per circle of radius 500m
simParams.NumRxMultiDU = simParams.NumRx*simParams.NumDU; % Number of receivers at the user end
simParams.posRx = [0;0;0];
% Configure Scatterers
simParams.refax = [[1;0;0] [0;1;0] [0;0;0]];

simParams.serveRadius = [20 30 40 50 100 200 500];
simParams.numServeRadius = length(simParams.serveRadius);
simParams.folderName = sprintf('MultiDUChannelModels_%dPaths',simParams.NumPaths);

%% === Configure the transmit and receive antenna elements for each pair of single DU
% and single user ===
simParams.txAntenna = phased.IsotropicAntennaElement;            % To avoid transmission beyond +/- 90
                                                                 % degrees from the broadside, baffle
                                                                 % the back of the transmit antenna
                                                                 % element by setting the BackBaffled
                                                                 % property to true.
                                                                 
simParams.rxAntenna = phased.IsotropicAntennaElement('BackBaffled',false); % To receive the signal from 360 degrees,
                                                                 % set the BackBaffled property to false

simParams.txArray = phased.NRRectangularPanelArray('Size',[1, 1, 1, 1],'ElementSet', {simParams.txAntenna},...
            'Spacing',[0.5*simParams.lambda,0.5*simParams.lambda,3*simParams.lambda,3*simParams.lambda]);
simParams.rxArray = phased.ULA('Element',simParams.rxAntenna, ...
    'NumElements',simParams.NumRx,'ElementSpacing',0.5*simParams.lambda,'ArrayAxis','x');
rxArrayStv = phased.SteeringVector('SensorArray',simParams.rxArray,'PropagationSpeed',simParams.c);

n = 0:simParams.NumRx-1;eta=pi;
% <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
% numFrames = 1; % 10 ms frames -> 1 subframe - 10 slots - 14 symbols
numUsefulOFDM = 12; % num subcarriers
numSCS = 12;
K = carrier.NSizeGrid * numSCS;
pilotLengths = [1 4 8 16 32 64];
numPilotLengths = length(pilotLengths);

%% Hold data separately for each DU in $$_SingleDU$$ data structures
NMSE_Hmk_MP_Pilot = zeros(numPilotLengths,1);
NMSE_Hmk_MP_NoPilot = zeros(numPilotLengths,1);
NMSE_Hmk_DFT_Pilot = zeros(numPilotLengths,1);
NMSE_Hmk_LMMSE_Pilot = zeros(numPilotLengths,1);

BER_MP_Pilot = zeros(numPilotLengths,1);
BER_MP_NoPilot = zeros(numPilotLengths,1);
BER_DFT = zeros(numPilotLengths,1);
BER_MMSE = zeros(numPilotLengths,1);
BER_True = zeros(numPilotLengths,1);

totalNumChannels = 30; % Monte-Carlo trials
totalNumTrialsPerChannel = 20;
tolFloating = 1e-2;
% 5.2129e-10 for 500 serving radius
Beta = 4.2129e-10; % Empirical Large-scale fading coefficient

%% >>>>>>>>>>>>>>>> MAIN SIMULATION LOOP -- Pilot Training <<<<<<<<<<<<<<<<<<<
snrdB = 10;snrPower = 10^(snrdB/10); % power
serveRadius = 100;
EbN0 = zeros(numPilotLengths,1);
simParams.totalNumSlots = 10;
simParams.totalNumChannels = totalNumChannels;
usefulOFDMRange = 1:numUsefulOFDM;
for channelIdx = 1:totalNumChannels
    for trialIdx = 1:totalNumTrialsPerChannel
        fprintf(">>>>>>>>>>>>>>>>>>> channelIdx=%d >>>>>>>>>>>>>>>>>>>>> \n", channelIdx)
        [AoAs_True,Amk_True,sigGridMultiDU,...
            noiseGridMultiDU,txGrid,txAmp] = genMultiDUChannelOutput(ofdmInfo,simParams,carrier,snrdB,serveRadius,channelIdx);
        sigGridMultiDU = sigGridMultiDU(:,usefulOFDMRange,:);
        noiseGridMultiDU = noiseGridMultiDU(:,usefulOFDMRange,:);
        rxGridMultiDU = sigGridMultiDU + noiseGridMultiDU;
        txGrid = txGrid(:,usefulOFDMRange);
        %% ========= Perform True Channel Estimation Algo. ===============
        Hest_Builtin = zeros(simParams.NumRxMultiDU,1);
        Hest_RandPhases = zeros(simParams.NumDU,1);
        for DUIdx = 1:simParams.NumDU
            antennaRange = (DUIdx-1)*8+1:DUIdx*8;
            sigGrid = sigGridMultiDU(:,:,antennaRange);
            noiseGrid = noiseGridMultiDU(:,:,antennaRange);
            rxGrid = sigGrid + noiseGrid;
            [hEstGrid,noiseEst] = nrChannelEstimate(rxGrid,txGrid);
            Hmk_true = reshape(hEstGrid(1,:,:),1,[],simParams.NumRx);
            Hmk_true = permute(Hmk_true,[3,2,1]); % NumRx X NumOFDMSymb
    
            ang_avg = 0;
            for i=1:12
                ang = angle(Hmk_true(1,i));
                ang_avg = ang_avg + ang;
            end
            ang_avg = ang_avg / 12;
            Hest_Builtin(antennaRange) = mean(Hmk_true,2);
            Hest_RandPhases(DUIdx) = ang_avg;
            % Hest_RandPhases(DUIdx) = angle( complex(Hest_Builtin((DUIdx-1)*8+1)) /  real(Hest_Builtin((DUIdx-1)*8+1)) );
        end
    
         %% ===========================================
         %  ===           Pilot Training            ===
         %  ===========================================
         pilotIdx = 1;
         while pilotIdx <= numPilotLengths
             
             currPilotL = pilotLengths(pilotIdx);
             freqMax = floor(currPilotL / numSCS) + 1;
             if mod(currPilotL,numSCS) == 1
                 currPilotL = 1;
             elseif freqMax == 1 && mod(currPilotL,numSCS) > 1
                 currPilotL = mod(currPilotL,numSCS);
             else
                 currPilotL = numSCS;
             end
    
             Hest_DFT = zeros(simParams.NumRxMultiDU,1);
             Hest_MP_Pilot = zeros(simParams.NumRxMultiDU,1);
             Hest_LinMMSE = zeros(simParams.NumRxMultiDU,1);
             Hest_MP_NoPilot = zeros(simParams.NumRxMultiDU,1);
    
             %% >>>>>>>>>>>>>> Channel Estimation at Single-DU Starts >>>>>>>>>>>>>>>>
             %txGrid:           NumSubcarrier X NumOFDMSym X
             %                       NumTransmitterAntenna
             %sigGridMultiDU:   NumSubcarrier X NumOFDMSym X
             %                                   (NumReceiverAntennaPerDU X
             %                                   NumDU)
             for DUIdx = 1:simParams.NumDU
                 antennaRange = (DUIdx-1)*8+1:DUIdx*8;
                 sigGrid = sigGridMultiDU(:,:,antennaRange);
                 noiseGrid = noiseGridMultiDU(:,:,antennaRange);
                 rxGrid = sigGrid + noiseGrid;
    
                 %% ========= Extract Signals from Resource Grid ===============
                 for freqIdx = 1:freqMax
                     Xpilot = txGrid(freqIdx,1:currPilotL);
                     YSampled = reshape(rxGrid(freqIdx,1:currPilotL,:),[],currPilotL,simParams.NumRx);     % pilot recv: NumOFDMSymb X NumRx
                     YSampled = permute(YSampled,[3,2,1]); % NumRx X NumOFDMSymb
                     NSampled = reshape(noiseGrid(freqIdx,1:currPilotL,:),[],currPilotL,simParams.NumRx);  % pilot recv: NumOFDMSymb X NumRx
                     NSampled = permute(NSampled,[3,2,1]); % NumRx X NumOFDMSymb
                     pilotNorm = norm(Xpilot);
    
                     %% ========= 1. Perform Angle-Domain Channel Estimation Algo. ===============
                     YTraining = YSampled*Xpilot'./norm(Xpilot)^2;
                     AoAs_dft = reshape(sort(dft_aoa(YTraining,simParams.NumRx,simParams.NumPaths),'ascend'),1,simParams.NumPaths);
                     AoAs_mp_pilot = sort(matpencil_aoa(YTraining,simParams.NumPaths),'ascend');
                     %  Channel Amplitude Estimation
                     Amk_dft = rxArrayStv(simParams.fc,[AoAs_dft;zeros(1,simParams.NumPaths)]);
                     Amk_mp = rxArrayStv(simParams.fc,[AoAs_mp_pilot;zeros(1,simParams.NumPaths)]);
                     Gmk_dft = Amk_dft*(Amk_dft\(YSampled*Xpilot'/pilotNorm^2));
                     Gmk_mp = Amk_mp*(Amk_mp\(YSampled*Xpilot'/pilotNorm^2));
    
                     Hest_DFT(antennaRange) = Hest_DFT(antennaRange) + Gmk_dft;
                     Hest_MP_Pilot(antennaRange) = Hest_MP_Pilot(antennaRange) + Gmk_mp;
    
                     %% ========= 2. Perform MMSE Channel Estimation  ===============
                     Hest_LinMMSE(antennaRange) = Hest_LinMMSE(antennaRange) + h_MMSE_CE(YSampled,Xpilot,NSampled,Beta);
   

                     tmp_mp_noPilot =  zeros(simParams.NumRx,1);
                     tmp_mp_noPilot2 = zeros(simParams.NumRx,1);
                     for ofdmSymIdx = 1:currPilotL
                         %% ========= 3. Perform Matrix-Pencil with No Pilot Training ===============
                         yAmp = norm(YSampled(:,ofdmSymIdx));
                         AoAs_mp_noPilot = sort(matpencil_aoa(YSampled(:,ofdmSymIdx),simParams.NumPaths),'ascend'); % zero-pilot
                         Amk_hat_mp_noPilot = rxArrayStv(simParams.fc,[AoAs_mp_noPilot;zeros(1,simParams.NumPaths)]);
                         tmp_mp_noPilot = Amk_hat_mp_noPilot*(Amk_hat_mp_noPilot\(YSampled(:,ofdmSymIdx)/(txAmp)));
    
                         % Remove (Random Phase + Random Symbol Phase) Altogether
                         phi_rand = angle(tmp_mp_noPilot(1));
                         tmp_mp_noPilot2 = tmp_mp_noPilot2 + tmp_mp_noPilot .* exp(-1i*phi_rand);
                     end
                     Hest_MP_NoPilot(antennaRange) = Hest_MP_NoPilot(antennaRange) + (tmp_mp_noPilot2 ./ currPilotL);
                 end
                 Hest_DFT(antennaRange) = Hest_DFT(antennaRange)./freqMax;
                 Hest_MP_Pilot(antennaRange) = Hest_MP_Pilot(antennaRange)./freqMax;
                 Hest_LinMMSE(antennaRange) = Hest_LinMMSE(antennaRange)./freqMax;
                 tmp = Hest_MP_NoPilot(antennaRange)./freqMax;
    
                 % angleOffset = angle( complex(Hest_MP_Pilot(antennaRange(1))) /  real(Hest_MP_Pilot(antennaRange(1))) );
                 angleOffset = Hest_RandPhases(DUIdx);
                 Hest_MP_NoPilot(antennaRange) = tmp .* exp(1i*angleOffset);
                 % Hest_MP_NoPilot(antennaRange) = tmp;
             end
             %% >>>>>>>>>>>>>> NMSE of Channel Estimation >>>>>>>>>>>>>>>
             nmseMP_Pilot = computeNMSE(Hest_Builtin,Hest_MP_Pilot);
             nmseMP_NoPilot = computeNMSE(Hest_Builtin,Hest_MP_NoPilot);
             nmseDFT = computeNMSE(Hest_Builtin,Hest_DFT);
             nmseLMMSE = computeNMSE(Hest_Builtin,Hest_LinMMSE);
    
             NMSE_Hmk_MP_Pilot(pilotIdx) = NMSE_Hmk_MP_Pilot(pilotIdx) + nmseMP_Pilot;
             NMSE_Hmk_MP_NoPilot(pilotIdx) = NMSE_Hmk_MP_NoPilot(pilotIdx) + nmseMP_NoPilot;
             NMSE_Hmk_DFT_Pilot(pilotIdx) = NMSE_Hmk_DFT_Pilot(pilotIdx) + nmseDFT;
             NMSE_Hmk_LMMSE_Pilot(pilotIdx) = NMSE_Hmk_LMMSE_Pilot(pilotIdx) + nmseLMMSE;
    
             YQPSKBER = reshape(rxGridMultiDU(1:K,usefulOFDMRange,:),[],K*numUsefulOFDM,simParams.NumRxMultiDU);     % pilot recv: NPilotSymb X NRx
             YQPSKBER = permute(YQPSKBER,[3,2,1]); % NumRx X NumQPSKSymb
             symEnc = txGrid(1:K,usefulOFDMRange) ./ txAmp;
             symEnc = reshape(symEnc,K*numUsefulOFDM,[]);
             berMP = computeBER(YQPSKBER,symEnc,Hest_MP_Pilot);
             berMP_NoPilot = computeBER(YQPSKBER,symEnc,Hest_MP_NoPilot);
             berDFT = computeBER(YQPSKBER,symEnc,Hest_DFT);
             berLMMSE = computeBER(YQPSKBER,symEnc,Hest_LinMMSE);
             berTrue = computeBER(YQPSKBER,symEnc,Hest_Builtin);
             BER_MP_Pilot(pilotIdx) = BER_MP_Pilot(pilotIdx) + berMP;
             BER_MP_NoPilot(pilotIdx) = BER_MP_NoPilot(pilotIdx) + berMP_NoPilot;
             BER_DFT(pilotIdx) = BER_DFT(pilotIdx) + berDFT;
             BER_MMSE(pilotIdx) = BER_MMSE(pilotIdx) + berLMMSE;
             BER_True(pilotIdx) = BER_True(pilotIdx) + berTrue;
    
             pilotIdx = pilotIdx + 1;
         end
     end
end
% NMSE_Hmk_MP_Pilot = NMSE_Hmk_MP_Pilot ./ totalNumChannels;
% NMSE_Hmk_MP_NoPilot = NMSE_Hmk_MP_NoPilot ./ totalNumChannels;
% NMSE_Hmk_DFT_Pilot = NMSE_Hmk_DFT_Pilot ./ totalNumChannels;
% NMSE_Hmk_LMMSE_Pilot = NMSE_Hmk_LMMSE_Pilot ./ totalNumChannels;
% 
% BER_MP_Pilot = BER_MP_Pilot ./ totalNumChannels;
% BER_MP_NoPilot = BER_MP_NoPilot ./ totalNumChannels;
% BER_DFT = BER_DFT ./ totalNumChannels;
% BER_MMSE = BER_MMSE ./ totalNumChannels;
% BER_True = BER_True ./ totalNumChannels;

figure
plot(unwrap(angle(Hest_Builtin)),'DisplayName','Theoretical');
hold on
plot(unwrap(angle(Hest_MP_NoPilot)),'DisplayName','No-Pilot Training');
legend show
grid on
hold off

NMSE_Hmk_MP_Pilot = NMSE_Hmk_MP_Pilot ./ (totalNumChannels*totalNumTrialsPerChannel);
NMSE_Hmk_MP_NoPilot = NMSE_Hmk_MP_NoPilot ./ (totalNumChannels*totalNumTrialsPerChannel);
NMSE_Hmk_DFT_Pilot = NMSE_Hmk_DFT_Pilot ./ (totalNumChannels*totalNumTrialsPerChannel);
NMSE_Hmk_LMMSE_Pilot = NMSE_Hmk_LMMSE_Pilot ./  (totalNumChannels*totalNumTrialsPerChannel);

BER_MP_Pilot = BER_MP_Pilot ./  (totalNumChannels*totalNumTrialsPerChannel);
BER_MP_NoPilot = BER_MP_NoPilot ./  (totalNumChannels*totalNumTrialsPerChannel);
BER_DFT = BER_DFT ./  (totalNumChannels*totalNumTrialsPerChannel);
BER_MMSE = BER_MMSE ./  (totalNumChannels*totalNumTrialsPerChannel);
BER_True = BER_True ./  (totalNumChannels*totalNumTrialsPerChannel);

%% >>>>>>>>>>>>>>>>>>>>> Plotting Start >>>>>>>>>>>>>>>>>>>>
job=string(datetime('now','Format',"yyyy-MM-dd-HH-mm-ss"));

fig1=figure;
semilogy(pilotLengths,NMSE_Hmk_MP_Pilot,'-o','DisplayName','Hmk Matrix Pencil (Pilot-Training)');
hold on;
semilogy(pilotLengths,NMSE_Hmk_MP_NoPilot,'-o','DisplayName','Hmk Matrix Pencil (Zero-Pilot)');
semilogy(pilotLengths,NMSE_Hmk_DFT_Pilot,'-^','DisplayName','Hmk DFT');
semilogy(pilotLengths,NMSE_Hmk_LMMSE_Pilot,'-*','DisplayName','Hmk MMSE');
grid on
xlabel('Number of Pilots')
ylabel('NMSE of CSI-Est')
title(sprintf('NMSE of CSI: NAnt=%d ServeRadius=%d NumPaths=%d',simParams.NumRx,serveRadius,simParams.NumPaths));
legend show
hold off
pngfile=sprintf('NMSE_CSI_zeroPilot_servR%d_%dPaths_j%s',serveRadius,simParams.NumPaths,job);
print(fig1,pngfile,'-dpng')


fig2=figure;
semilogy(pilotLengths,BER_MP_Pilot,'-o','DisplayName','BER Matrix-Pencil (Pilot-Training)');
hold on
semilogy(pilotLengths,BER_MP_NoPilot,'-o','DisplayName','BER Matrix-Pencil (Zero-Pilot)');
semilogy(pilotLengths,BER_DFT,'-^','DisplayName','BER DFT');
semilogy(pilotLengths,BER_MMSE,'-*','DisplayName','BER Linear-MMSE');
grid on
xlabel('Number of Pilots')
ylabel('BER')
title(sprintf('BER vs Number of Pilots ServeRadius=%d NumPaths=%d',serveRadius,simParams.NumPaths));
legend show     
hold off
pngfile=sprintf('BER_CSI_zeroPilot_servR%d_%dPaths_j%s',serveRadius,simParams.NumPaths,job);
print(fig2,pngfile,'-dpng')


% fig3=figure;
% hold on
% semilogy(pilotLengths,EbN0,'-o','DisplayName','EbN0');
% grid on
% xlabel('Number of Pilots')
% ylabel('EbN0')
% title(sprintf('Num. of Pilots vs EbN0 NumDU=%d ServeRadius=%d NumPaths=%d',simParams.NumDU,serveRadius,simParams.NumPaths));
% legend show
% hold off
% pngfile=sprintf('EbN0_Pilot_MultiDU_servR%d_%dPaths_j%s',serveRadius,simParams.NumPaths,job);
% print(fig3,pngfile,'-dpng')
%% <<<<<<<<<<<<<<<<<< Plotting End <<<<<<<<<<<<<<<<<<<<<<<<
%% =============== END OF MAIN FUNCTION ==================


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                   HELPER FUNCTIONS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [AoAs_True_MultiDU,Amk_True_MultiDU,sigGridMultiDU,noiseGridMultiDU,txGrid,txAmp] = genMultiDUChannelOutput(...
    ofdmInfo,simParams,carrier,snrdB,serveRadius,channelIdx)
refax = [[1;0;0] [0;1;0] [0;0;0]];
% Transmitter Setup
txAmp = 3000*10^(snrdB/20); % 40dBm
nSlot = 0;
carrier.NSlot = nSlot;
slotGrid = nrResourceGrid(carrier,simParams.NumTx);
qamSymbols = nrSymbolModulate(randi([0,1],numel(slotGrid)*2,1),'QPSK');
slotGrid(:) = txAmp.*qamSymbols;
txGrid= slotGrid; % 12 x 14

% Multi-DU Setup
sigGridMultiDU = [];
AoAs_True_MultiDU = [];
Amk_True_MultiDU = [];
rxArrayStv = phased.SteeringVector('SensorArray',simParams.rxArray,'PropagationSpeed',simParams.c);
for DUIdx = 1:simParams.NumDU
    chanFileName = fullfile(simParams.folderName,sprintf('R%d-Chan%d-DU%d.mat', ...
                    serveRadius,channelIdx,DUIdx));
    file = java.io.File(chanFileName);
    fullpath = char(file.getAbsolutePath());

    chanModel = load(fullpath);
    
    channel = chanModel.channel;
    simParams.scatPos = channel.ScattererPosition;
    simParams.posTx = channel.TransmitArrayPosition;
    simParams.posRx = channel.ReceiveArrayPosition;
    [~,~,tau] = channel(complex(randn(ofdmInfo.SampleRate*1e-3,simParams.NumTx), ...
        randn(ofdmInfo.SampleRate*1e-3,simParams.NumTx)));
    maxChDelay = ceil(max(tau)*ofdmInfo.SampleRate);
    
    
    %% === Send the Waveform through the Channel ===
    % OFDM Modulation
    [txWaveform,~] = nrOFDMModulate(carrier,txGrid);
    % Append zeros to the transmit waveform to account for channel delay
    txWaveform = [txWaveform; zeros(maxChDelay,simParams.NumTx)];
    % Pass the waveform through the channel
    [fadWave,~,~] = channel(txWaveform);
    
    % Estimate timing offset
    %% >>>>>>>>>>>>>>>>>> Channel Estimation Start >>>>>>>>>>>>>>>>>>>>
    offset = nrTimingEstimate(carrier,fadWave,txGrid);
    if offset > maxChDelay
        offset = 0;
    end
    
    % Receiver Setup
    
    % Compute True AoA based on ScatPos and posRx
    [~,AoAs_True] = rangeangle(simParams.scatPos,simParams.posRx,refax);
    AoAs_True(1,:) = sort(AoAs_True(1,:),'ascend');
    Amk_True = rxArrayStv(simParams.fc,AoAs_True(1,:));
    
    
    % Correct timing offset
    fadWave = fadWave(1+offset:end,:);
    % Perform OFDM demodulation
    sigGrid = nrOFDMDemodulate(carrier,fadWave);           
   
    AoAs_True_MultiDU = cat(2,AoAs_True_MultiDU,AoAs_True);
    Amk_True_MultiDU = cat(2,Amk_True_MultiDU,Amk_True);
    sigGridMultiDU = cat(3,sigGridMultiDU,sigGrid);
end
noiseGridMultiDU = 1/sqrt(2).*complex(randn(size(sigGridMultiDU)),randn(size(sigGridMultiDU)));
end


%% >>>>>>>>>>>>>>>  NMSE >>>>>>>>>>>>>>>
function nmse = computeNMSE(H_true,H_est)
% Input:
%   H_true / H_est: NumRx x NumOFDMSym (Time Domain Samples)
    sqErr = mean(abs(H_true - H_est).^2);
    sqNorm_H_true = mean(abs(H_true).^2);
    nmse = sqErr / sqNorm_H_true;
end
%% >>>>>>>>>> BER >>>>>>>>>>>>
function BER = computeBER(yQPSK,symEnc,Hmk_est)
    [~,nSymb] = size(yQPSK);
    symDec = Hmk_est'*yQPSK;
    symEncQPSK = nrSymbolDemodulate(symEnc,'QPSK','DecisionType','Hard');
    symDecQPSK = nrSymbolDemodulate(symDec.','QPSK','DecisionType','Hard');
    BER = biterr(symEncQPSK,symDecQPSK)/(2*nSymb);
end

%%======== Channel Estimation ========
function h_LinMMSE = h_MMSE_CE(y,x,n,Beta)
% y                 = Frequency-domain received signal
%                   NumRx X NumSym
% d                 = Large-scale fading parameter
%                   NumRx X NumPaths
% x                 = pilot symbol
% noise
pilotNorm = norm(x);
[NumRx,currPilotL] = size(y);
% h_LinMMSE = zeros(NumRx,1);
% for p=1:currPilotL
%     x_sqval = sum(abs(x(p)).^2);
%     yExtracted = y(:,p)*x(p)'/x_sqval;
%     zeta = Beta / Beta + (trace(n(:,p)*n(:,p)')/trace(y(:,p)*y(:,p)'));
%     W_MMSE = 1/(1+1/zeta);
%     h_LinMMSE = h_LinMMSE + W_MMSE*yExtracted;
% end
% h_LinMMSE = h_LinMMSE ./ currPilotL;
yExtracted = y*x'/(pilotNorm^2);
zeta = Beta / Beta + (trace(n*n')/trace(y*y')/currPilotL);
W_MMSE = 1/(1+1/zeta);
h_LinMMSE = W_MMSE*yExtracted;
end



%% AoA Estimation Related
%% ===== Extended DFT + Angle-Rotation =====
function [AoA_DFT] = dft_aoa(ymk_Sampled,N,L)
% Input:    yt, NumRx, N_MultiPath
% Output:   AOA_estimated, beta_estimated
[~,T] = size(ymk_Sampled);
FN = dftmtx(N)/sqrt(N);
Ndft_points = 100; %% can choose whatever # you want
AoA_DFT = zeros(L,1);

for t=1:T
    AoA = zeros(L,1);
    hDFT = FN * ymk_Sampled;
    
    % Coarse Peak Finding
    % -- Find the central point (qInit) of each path
    [qInits,isNeg] = findInitDFTBin(hDFT,N,L);
    
    for l=1:L
        fNq = FN(:,qInits(l));
        ymk_DFT = fNq .* ymk_Sampled;
        
        angles_in_phi = [-Ndft_points/2: Ndft_points /2]*pi/ Ndft_points; %% Ndft_points in the phi domain
        st_vec_mtx = exp(1i* [0:N-1]' * angles_in_phi);  %% N \times Ndft_points matrix of Ndft_points steering vectors
    
        % Now if x is the data vector
        angle_info = abs(st_vec_mtx' * ymk_DFT);
        [~, max_angle_location] = max(angle_info);
        phi_location = angles_in_phi(max_angle_location);
        theta_init = 2*qInits(l)/N;
        
        if isNeg(l)
            theta = -theta_init + phi_location/pi; %% since \phi = kd\sin\theta
        else
            theta = theta_init - phi_location/pi;
        end
    
        if abs(theta) > 1
            theta = findNextPeak(theta,qInits(l),angle_info,angles_in_phi,N);
        end
    
        if isNeg(l)
            AoA(l) = -1*real(acosd(theta));
        else
            AoA(l) = real(acosd(theta));
        end
    end
    AoA = sort(AoA);
    AoA_DFT = AoA_DFT + AoA;
end
AoA_DFT = AoA_DFT ./ T;
end


function [Q,isNeg] = findInitDFTBin(hDFT,N,L)
    [~,I] = sort(abs(hDFT),'descend');
    threshold = floor(N/2);
    Q = zeros(1,L);
    isNeg = zeros(1,L);
    pl = 1;
    
    for l=1:N
        if I(l) >= (threshold + 1)
            Q(pl) = I(l)-2;
            pl = pl+1;
            isNeg(l) = 1;
        else
            Q(pl) = I(l);
            pl = pl+1;
            isNeg(l) = 0;
        end
        if pl > L
            break;
        end
    end
end

function new_theta = findNextPeak(prev_theta,curr_qInit,angle_info,angles_in_phi,N)
    [ang_pks,ang_loc] = findpeaks(angle_info);
    [~,sorted_ang_loc_ind] = sort(ang_pks,'descend');
    ang_loc_sorted = ang_loc(sorted_ang_loc_ind);
    
    ang_locs_L = ang_loc_sorted(2:end);
    isNeg = 0; 
    new_theta = 2; % dummy init value
    if sign(prev_theta) < 1
        isNeg = 1;
    end
    idx=1; NPks = length(ang_locs_L);
    while (abs(new_theta) > 1 && idx <= NPks)
        curr_max_angle_loc = ang_locs_L(idx);
        curr_phi_loc = angles_in_phi(curr_max_angle_loc);
        if isNeg
            new_theta = -2*curr_qInit/N + curr_phi_loc/pi;
        else
            new_theta = 2*curr_qInit/N - curr_phi_loc/pi;
        end
        idx = idx+1;
    end
end


%% ====== AoA estimation using Matrix Pencil ======
function AoA_MP = matpencil_aoa(ymk_Sampled,L)
% ymk_Sampled = NxT
% y = As + n
% N: Number of array elements
% M: Number of paths
% L: Matrix Pencil parameter
% T: Number of time samples

[N,T] = size(ymk_Sampled);
AoA_MP = zeros(1,L);
for t=1:T
    P = L; % size of window
    x = ymk_Sampled(:,1); % Kx1=>only one time sample
    Y1 = zeros(N-P,P);
    Y2 = zeros(N-P,P);
    for p=1:P
        Y1(:,p) = x(p:N-P-1+p,1);
        Y2(:,p) = x(p+1:N-P+p,1);
    end
    
    Y1_pinv = (Y1'*Y1)\Y1';
    z_hat = sort(eig(Y1_pinv*Y2),'descend');
    
    % AOAs = sort(acosd(imag(log(z_hat(1:L))) ./ pi));
    % AOAs = AOAs';
    
    % find the upper bound (need |e^i pi cos(theta)| <=1 )
    [~,i] = mink(abs(abs(z_hat)-1),L);
    z_hat_cap = z_hat(i);
    AOAs = zeros(1,L);
    
    for pl=1:L
        z_tmp = z_hat_cap(pl);
    
        theta_tmp = acosd(imag(log(z_tmp)) ./ pi);
        theta_tmpReal = asind(imag(log(z_tmp)) ./ pi);
        if sign(theta_tmpReal) > 0
            AOAs(pl) = theta_tmp;
        else
            AOAs(pl) = -theta_tmp;
        end
    end
    AoA_MP = AoA_MP + sort(AOAs);
end
AoA_MP = AoA_MP ./ T;
end
