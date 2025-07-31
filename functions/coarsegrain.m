%%% coarsegrain
%
% Perform temporal coarsegraining of univariate or multichannel time series.
%
% This function downsamples a signal by averaging over fixed-size time windows. 
% It supports both vector and matrix input (channels × time), and can operate 
% using either non-overlapping or overlapping windows depending on the step size.
%
% Units for the window and step can be specified in data points, milliseconds, 
% or seconds, with optional sampling rate conversion.
%
% --- INPUTS ---
%   TS           : [C × T] matrix or [1 × T] vector of time series data  
%                  If matrix, rows are assumed to be channels and columns time points.
%
%   Win          : Window size for coarsegraining (in points, seconds, or milliseconds)
%
% --- Optional Parameters (name-value pairs) ---
%   'step'       : Step size between successive windows (default = Win)
%
%   'unit'       : Unit used for Win and step. Options:
%                    - 'pts' (default) : number of data points
%                    - 's'             : seconds
%                    - 'ms'            : milliseconds
%
%   'Fs'         : Sampling frequency (in Hz). Required if unit is 's' or 'ms'
%   'coarsetype' : Calculus used to coarsegrain, can be 'mean' or 'median'. Default is 'mean'
%
% --- OUTPUTS ---
%   Coa_Sig      : Coarsegrained signal matrix [C × N]
%                  Each row corresponds to a channel, each column to a window.
%
% --- Notes ---
% - Input matrix must be formatted as [channels × time]. Use transpose if necessary.
% - If 'unit' is time-based ('s' or 'ms'), the sampling rate must be provided via 'Fs'.
% - If TS is a row or column vector, it is automatically reshaped to 1 × T.
% - Output matrix Coa_Sig will have as many columns as there are full windows.
%
% --- Example Usage ---
%   fs = 1000;                                 % Sampling rate (Hz)
%   ts = randn(4, 10000);                      % 4-channel signal, 10 s long
%   coarse = coarsegrain(ts, 200, ...
%                        'step', 100, ...
%                        'unit', 'pts', ...
%                        'Fs', fs);            % 200-pt window, 100-pt step
%
%   coarse_t = coarsegrain(ts, 0.2, ...
%                          'step', 0.1, ...
%                          'unit', 's', ...
%                          'Fs', fs);          % 200 ms window, 100 ms step
%
% --- Author ---
% Matthieu Aguilera — FunSy team, LNCA (Strasbourg), July 2025

function Coa_Sig = coarsegrain(TS, Win, varargin)

% Parse inputs
p = inputParser;
addParameter(p, 'step', nan);
addParameter(p, 'unit', 'pts');
addParameter(p, 'Fs', nan);
addParameter(p, 'coarsetype', 'mean');
parse(p, varargin{:});
In = p.Results;

% Convert window to points
if ~strcmp(In.unit, 'pts') && isnan(In.Fs)
    error('Please provide Sampling Frequency Fs if you use "s" or "ms" as unit.');
end

switch In.unit
    case 's'
        W = round(Win * In.Fs);
    case 'ms'
        W = round((Win / 1000) * In.Fs);
    case 'pts'
        W = Win;
    otherwise
        error('Unknown unit type. Use ''pts'', ''s'' or ''ms''.');
end

% Validate dimensions
if size(TS,1) > 1 && size(TS,1) > size(TS,2)
    error('Input must be in channels x time format.');
end

% Ensure TS is 2D: row vector becomes 1×N
if isvector(TS)
    TS = reshape(TS, 1, []);  % Row vector
end

L = size(TS,2);
if isnan(In.step)
    step = W;
else
    step = In.step;
end

if L < W
    error('Window size is too large for the input signal length.');
end

q = floor((L - W) / step) + 1;
nChannels = size(TS,1);
Coa_Sig = zeros(nChannels, q);

for i = 1:q
    idx_start = (i-1)*step + 1;
    idx_end = idx_start + W - 1;
    switch In.coarsetype
        case 'mean'
        Coa_Sig(:,i) = mean(TS(:,idx_start:idx_end), 2);
        case 'median'
        Coa_Sig(:,i) = median(TS(:,idx_start:idx_end), 2);
    end
end
