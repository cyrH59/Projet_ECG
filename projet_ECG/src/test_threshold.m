%% Main script to test ecg function without gui
% This file computes a simple analysis of an ecg signal. You can use it to test the different processing methods. 
% This first version will plot the temporal signal, compute its cardiac rythma and display the different P, Q, R, S, T points for a specific segment.  

clear; close all; clc;
addpath(genpath('.'));

%% Load a signal
%[file,path] = uigetfile('*.mat', 'rt');
%signal = load(fullfile(path, file));
signal=load("../data/ecg_normal_1.mat");
signal2=load("../data/ecg_VF.mat");
data = -signal.ecg;
data2=-signal2.ecg 
% Your ecg data
Fs = signal.Fs; % Sampling frequency
Fs2=signal2.Fs
N = size(data,2); % Data length
time_axis = (1:N)/Fs;

%% Threshold method
th = 200; % threshold
i_seg = 10; % Segment number to plot

% Time plot
figure;
plot(time_axis, data); grid on;
hold on; plot(time_axis, th*ones(1,N), 'red');
xlabel('Time (s)');
ylabel('Magnitude');
title('Time evolution of the loaded signal')

% Print BPM
[bpm, R_locs] = bpm_threshold(data, th, Fs);
% Figures PQRST
[segment, P_loc, Q_loc, R_loc, S_loc, T_loc] = ecg_threshold(data, R_locs, i_seg);
time_segment = (1:length(segment))/Fs;

figure;
h = plot(time_segment, segment); grid on;
hold on;
plot(time_segment(P_loc),segment(P_loc), '*','Color','red'); text(time_segment(P_loc),segment(P_loc),' P ','Color','red','FontSize',14);
plot(time_segment(Q_loc),segment(Q_loc), '*','Color','red'); text(time_segment(Q_loc),segment(Q_loc),' Q ','Color','red','FontSize',14);
plot(time_segment(R_loc),segment(R_loc), '*','Color','red'); text(time_segment(R_loc),segment(R_loc),' R ','Color','red','FontSize',14);
plot(time_segment(S_loc),segment(S_loc), '*','Color','red'); text(time_segment(S_loc),segment(S_loc),' S ','Color','red','FontSize',14);
plot(time_segment(T_loc),segment(T_loc), '*','Color','red'); text(time_segment(T_loc),segment(T_loc),' T ','Color','red','FontSize',14);
hold off;
xlabel('Time (s)');
ylabel('Magnitude');
title('ECG segment characteristic')

%% Your turn : My new method !

figure()
x=data;
w=ones(8*Fs,1);
d=300;
N_fft=100;
[Sx, f, t] = spectro(x,w,d,N_fft,Fs);
imagesc(f,t,20*log(Sx));
title("Spectrogramme ecg normal 1");
xlabel("Fréquence (Hz)");
ylabel("Temps (s)");

figure()
x2=data2;
w2=ones(8*Fs2,1);
d2=300;
N_fft2=100;
[Sx2, f2, t2] = spectro(x2,w2,d2,N_fft2,Fs2);
imagesc(f2,t2,20*log(Sx2));
title("Spectrogramme ecg VF");
xlabel("Fréquence (Hz)");
ylabel("Temps (s)");

