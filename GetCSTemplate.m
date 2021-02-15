function [psi, TECGShift] = GetCSTemplate(peaks,ecg,Fs,N,m,Phi,EstType)
%
% [psi, TECGShift] = GetCSTemplate(refQRS,ecg,Fs,WinLength,NumMeasu,Phi,EstType)
%
% Overview: return a matrix containing in each row a a compressed translated 
%           version of the 'reference' QRS template and the uncompressed
%           version
% Inputs:
%       peaks : location of QRSs to be used for building the QRS template
%         ecg : original ecg signal
%          Fs : sampling frequency
%           N : Number of sample of a segment before compression
%           m : number of measurements for compressive sensing
%         Phi : sensing matrix
%     EstType : type of compressed sensed estimator 'direct' or 'orth'
%
% Outputs:
%         psi   : matrix of compressed shifthed template
%     TECGShift : matrix of shifted template
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

% ---  Generate Template for QRS ----

TECG = BuildTemplate(peaks, ecg, Fs);    % Build the average template

% Shift the template for compression pourposes 
TECGpad=[ TECG zeros(1,N)];                 % Zero padding to match the length of the window used 
TECGShift=gallery('circul',TECGpad);        % Generate a matrix with all the possible shift within the vector length
TECGShift=TECGShift(1:N,19:19+N-1);  

% Remove the mean from each shifted template 
for i=1:N
    TECGShift(i,:)=TECGShift(i,:)-mean(TECGShift(i,:));
end


% ---  Compressed Estimator Selection and Computation ---

switch EstType
    case 'direct'
    % Direct estimator
    psi = (Phi*TECGShift')';  % each row of psi is a compressed translated version of the template
    case 'orth'
    % Orthogonalized estimator
    InvOp = inv(Phi*Phi')*Phi;
    psi = (N/m*InvOp*TECGShift')';  % each row of psi is a compressed translated version of the template
end
