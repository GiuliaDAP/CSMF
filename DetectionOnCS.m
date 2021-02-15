function QRS = DetectionOnCS(y, psi, k, N, Fs, Phi)
%
% QRS = DetectionOnCS(y, psi, k, N, Fs , Phi)
%
% Overview
%    Detect QRS complexes (R peaks) on a compressed sensed signal
%   
% Inputs:      
%       y    : matrix (kxN), compressed sensed signal, each raw is a 
%              compressed segment of ECG 
%       psi  : matrix (NxN), compressed version of the shifted template, 
%       k    : number of CS blocks to be analyzed for R-peak detection                   
%       N    : Length of each original block (ECG segment)
%       Fs   : Sampling frequency
%       Phi  : mxN matrix, Sensing matrix used to compress the signal
% Outputs:
%       QRS  : location of detected QRS complexes 
%
% Reference: 
%   G Da Poian, CJ Rozell, R. Bernardini, R Rinaldo and GD Clifford, 
%   "Matched Filtering for Heart Rate Estimation on Compressive Sensing
%   ECG Measurements," in IEEE Transactions on Biomedical Engineering, 2017
%   doi: 10.1109/TBME.2017.2752422
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


thL = 0.5;     % First threshold level 
thH = 0.75;    % Second threshold level



QRS = [];
CorrVal = zeros(k,N); 
th = zeros(k,1);
val = zeros(k,1);
pos = cell(k);
nobeats = zeros(k,1);
idx = 1;
InvOp = inv(Phi*Phi')*Phi;
m=size(Phi,1);

% Detect QRS on each compressed segment
while idx<=k
    
    % Mean value estimation from compressed measurements 
    EstMean = y(idx,:)*(N/m*InvOp*1/N*ones(N,1));
   
    % Remove the estimated mean from the compressed signal
    y(idx,:) = y(idx,:)-(Phi*(EstMean*ones(N,1)))'; 
       
    % Correlation Estimation
    CorrVal(idx,:) =  abs(y(idx,:)*psi');  
    
    % Compute the adaptive threshold level
    [val(idx),~] = max(CorrVal(idx,:));
    
    if (rms(CorrVal(idx,:))) > 0.25*val(idx)  
       th(idx) = thH*val(idx); 
    else
       th(idx) = thL*val(idx);    
    end
    
   
    if idx>1 && val(idx)<1/4*val(idx-1) 
       pos{idx} = [];
       nobeats(idx) = 1;
    else   
       % look for three peaks instead of one
       [~,pos{idx}] = findpeaks(CorrVal(idx,:),'MinPeakDistance',...
                      round(Fs*0.25), 'MinPeakHeight',th(idx));    
       nobeats(idx)=0;
    end
    
    pos{idx} = pos{idx}+(idx-1)*N;
    
   % look for double peaks between to consecutive block and missing peak

   if ~isempty(QRS) && idx<k-1
      actualQRS = cell2mat(pos(idx));
      if ~isempty(actualQRS)
           if (abs(QRS(end)-actualQRS(1))< 0.2*Fs || abs(QRS(end)-actualQRS(1))< 0.2*medianRR )  % double peak
               if length(actualQRS) >1          
                    QRS = [QRS(1:end-1) round(mean([ QRS(end) actualQRS(1)])) actualQRS(2:end)];
               elseif length(actualQRS) == 1 
                    QRS = [QRS(1:end-1) round(mean([ QRS(end) actualQRS(1)])) ];
               end
           elseif abs(QRS(end)-actualQRS(1)) > 1.8*medianRR && idx>3  %  missed beat between two consecutive blocks
               [~,missPeak] = findpeaks([CorrVal(idx-1,N-40:N) ...
                              CorrVal(idx,1:40)], 'MinPeakHeight',th(idx)*0.75,'NPeaks',1); 
                              % look around a window of +/- 40 samples on
                              % the border and lower the th of 25%
                QRS = [QRS missPeak+N-40+(N*(idx-2)) actualQRS];
           else
             QRS = [QRS actualQRS];
          end 
      end
   else
       QRS = [QRS cell2mat(pos(idx))];
   end
   
   if idx <= 10
      medianRR = median(diff(QRS));
   else
      % QRS in last 10 blocks
      tmpQRS = QRS( QRS> N * (idx-10));
      medianRR = median(diff(tmpQRS));
   end

   idx=idx+1;

end



