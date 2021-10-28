function [X, f, t] = stft(x,w,d,N_fft,Fs)
% This function computes the stft for m = [0, d, 2d, 3d...]
% This function outputs are:
% -> X, which is  a matrix of n_fft lines and M columns
%    M is the number of elements of m
%    X(i,j) is the value of the spectrogram for time t(i) and frequency f(j)
% -> f, is a column vector of the frequencies (in Hz)
% -> t, is a row vector containing the times of the beginning of the windows
N = 1600 ;
m = zeros(N,fix(length(x)/d)+1);
for i = 1 : fix(length(x)/d)+1
    m(1:min(N,length(x)-(d*(i-1))),i) =  x((i-1)*d+1:min(length(x),(i-1)*d+N))';
    m(:,i) = m(:,i).*w ;
end
Ts = 1/Fs;
t = (0:d*Ts:Ts*length(x));
X = fft(m,N_fft);
f = (0:1/N_fft*Fs:(N_fft-1)/N_fft*Fs)';