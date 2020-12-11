function [Y,f] = easyFFT( X, n, dim, samplerate )
%easyFFT Easy-to-use Fast Fourier Transform
%  Conveniently returns the frequency vector along with the FFT. 
%
%  Input Arguments
%   - SIG        : Input array, signal in time domain, equidistant sampling 
%   - n          : Transform length 
%   - dim        : Dimension to operate along
%   - samplerate : Samplig rate of X in Hz 
%
%  Output Arguments
%   - Y          : Frequency domain representation (fftshift'ed format) 
%   - f          : Frequency vector along dim
%
%  Usage Example
%   >> [Y,f] = easyFFT( X, 2^20, 1, 1/dt );
%
% arguments  
%   X                double {mustBeNumeric,mustBeReal}        % Input array
%   n          (1,1) double {mustBeInteger,mustBeNonnegative} % Transformation length TODO accept also []  
%   dim        (1,1) double {mustBeInteger,mustBePositive   } % Dimension to operate along
%   samplerate (1,1) double {              mustBePositive   } % sampling rate of X 
% end

% Perform FFT
% Y = fft(X,n,dim); 
% Y = fftshift(Y,dim); 
Y = fft(X); 
Y = fftshift(Y); 
% Determine frequency vector
nSamples = length(Y); % number of samples
df = samplerate / nSamples;
f = 0 : df : df*(nSamples-1);
f = fftshift(f);
ind = f>=samplerate/2-eps(samplerate/2);
f(ind) = f(ind)-samplerate;
f = shiftdim(f, 2-dim); 

end


