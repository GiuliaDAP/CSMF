function TECG = BuildTemplate(peaks, ecg, Fs)

% TECG = BuildTemplate(peaks, ecg, Fs)
%
% Overview: QRS template computed by taking the median of the QRS complexes
%           in the ecg signal 
%
% Inputs:      
%       peaks  : Location of the R-peaks
%       ecg    : Original ecg signal
%       Fs     : Sampling frequency
% Outputs:
%       TECG   : QRS template 
%      
% Authors
%    Giulia Da Poian <giulia.dap@gmail.com>
%
% Reference: 
%   G Da Poian, CJ Rozell, R. Bernardini, R Rinaldo and GD Clifford, 
%   "Matched Filtering for Heart Rate Estimation on Compressive Sensing
%   ECG Measurements," in IEEE Transactions on Biomedical Engineering, 2017
%   doi: 10.1109/TBME.2017.2752422
%
% Copyright (C) Authors, all rights reserved.
%
% This software may be modified and distributed under the terms
% of the BSD license.  See the LICENSE file in this repo for details.
% This software may be modified & distributed under the terms
% of the BSD license. See LICENSE file in repo for details.
% Isolate days in this data

Pstart = round(0.05*Fs); 
Tstop = round(0.05*Fs);  
 
NB_CYCLES = length(peaks)-2;
% == template ecg (TECG)
indMQRSpeaks = find(peaks>Pstart);


for  cc=1:NB_CYCLES
    peak_nb = peaks(indMQRSpeaks(cc+1));   % +1 to unsure full cycles
    ecg_buff(:,cc) = ecg(peak_nb-Pstart:peak_nb+Tstop)';
end
TECG = median(ecg_buff,2);

TECG = (TECG-mean(TECG))';

