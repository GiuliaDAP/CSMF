% Demo for matched filtering approach in the compressd sensed domain
%
% Reference:
%         "G Da Poian, CJ Rozell, R. Bernardini, R Rinaldo and GD Clifford, 
%          "Matched Filtering for Heart Rate Estimation on Compressive 
%          Sensing ECG Measurements," in IEEE Transactions on Biomedical 
%          Engineering, 2017 doi: 10.1109/TBME.2017.2752422"
%
% Authors
%    Giulia Da Poian <giulia.dap@gmail.com>
%
% Copyright (C) Authors, all rights reserved.
%
% This software may be modified and distributed under the terms
% of the BSD license.  See the LICENSE file in this repo for details.
% This software may be modified & distributed under the terms
% of the BSD license. See LICENSE file in repo for details.
% Isolate days in this data

clear
clc
warning('off','all') 

Fs = 360;  % Sampling frequency

% Setting parameters for compression
WinLength  = 540;  % Window length (should be ~1.5 s of signal)
NumMeasu = 540/4;    % CR = 75% NumMeasu = 128

% Compressed Detector settings 'orth' or 'direct'
EstType = 'orth';

% Traditional Detector settings
uncomp_qrs_tresh = 0.6;
sec_temp = 10; % Build template on 10 sec of ecg signal


%%% START %%%

% 1. Load signal 
fid=fopen('test_ecg.dat','r');
f=fread(fid,2*Fs*10000,'bit12');
ecg=f(1:2:length(f));

Nwin = floor(length(ecg)/WinLength);  % Compute the maximum number of windows to be analyzed

% 2. Signal sgementation in non overlapping windows of length WinLength ~ 1.5s each 
ecg = ecg(1:WinLength*floor(length(ecg)/WinLength));
SegmSignal = reshape(ecg,WinLength,floor(length(ecg)/WinLength)); % SegmSignal is a [WinLengthxNwin] signal matrix

% 3. Signal compression using Compressive Sensing 
Phi = randn(NumMeasu,WinLength)./sqrt(NumMeasu);  % Create the sensing matrix 
CompSignalMatrix = (Phi*SegmSignal)';  


% 4. Generate a compressed template
% --- Need to get uncompressed QRS detection for initialization 
uncomp_QRS = qrs_detect2(ecg(1:Fs*sec_temp),0.250,uncomp_qrs_tresh,Fs,[],0);
% --- use only first XX s of signal for template generation
QRSforTemp = uncomp_QRS(uncomp_QRS <= Fs*sec_temp);
% --- generate compressed template
[psi, TECG] = GetCSTemplate(QRSforTemp,ecg(1:sec_temp*Fs), ...
                                    Fs,WinLength,NumMeasu,Phi,EstType);

% 5. R wave detection on compressed measurements ---
QRS_onCS = DetectionOnCS(CompSignalMatrix, ...
               psi,Nwin,WinLength, Fs,Phi);  % Compressed peak detector    

           
plot(ecg)
hold on
plot(QRS_onCS, ecg(QRS_onCS),'x')
legend('ECG', ['QRS detected using CSMF @ CR = ' num2str(round(100*(1-(NumMeasu/WinLength)))) '%'])
grid on
                                