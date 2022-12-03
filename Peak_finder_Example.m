%% proposed algorithm for detecting peaks in spectrum of brian signal
%% Ashkan Farrokhi
%% Date: '29-Nov-2022'
%% IUST
%% 
clear
clc
close all hidden
%% Load Signal:
load 'LFP_Simulated_Signal.mat'
S1 = Report.Signal.S1;                                                     % Source 1
S2 = Report.Signal.S2;                                                     % Source 2
S3 = Report.Signal.S3;                                                     % Source 3
Signal = Report.Signal.Observation;                                        % Sum(Sources)+Background
Fs=Report.Propety.Fs;

%=based on Our paper:
NumberOfComponent=3;                                                       %How many peaks (NOP)/ output of function 
MaxPointer=7;                                                              %Maximum number of peak extracted during calculation 
SpecOption.L_Segment=3;                                                    %(Second), window  size used by IRASA >>> when calculating fractal compunent
SpecOption.Overlap=.4;                                                     %(Second), Overlap size used by IRASA >>> when calculating fractal compunent
SpecOption.Freq_Lim=[4 200];                                               %(Hz), Frequency Limitation for calculation

%= Fvec: frequency range array 
%= frequency ranges that calculated frequency peaks should be between them
%= two way for define this parameter:
%= first: defined a [1*2] array >> Fvec(1): start frequency ,  Fvec(2): stop frequency
%= example: Fvec=[4 50]; >>> the algorithm only accept peaks that would be
%between 4 and 50 Hz
Fvec=[4 50];
%= first: defined a [NOP*2] array >> Fvec(NOP,1): start frequency for NOPth peak ,  Fvec(NOP,2): stop frequency for NOPth peak
%= example: Fvec=[8 12; 18 22; 28 32]; >>> the algorithm search for one
%peak inside of each of ranges >> [8 12], [18 22], [28 32]
% % Fvec=[8 12; 18 22; 28 32];

%=calculate peaks >> f0Vev (NOP,1)
[f0Vec, Width]= AFsearchForF0(Signal,Fs,NumberOfComponent,MaxPointer,SpecOption,Fvec);
disp(f0Vec)

%% plot: ... PSD & Location of peaks 
figure(1)


figure(2), hold on
title('PSD')
[pxxS,f]=pwelch(Signal,[],[],[],Fs);
plot(f,10*log10(pxxS),'k','LineWidth',3)

[pxx,f]=pwelch(S1,[],[],[],Fs);
plot(f,10*log10(pxx),'b','LineWidth',2)

[pxx,f]=pwelch(S2,[],[],[],Fs);
plot(f,10*log10(pxx),'b','LineWidth',2)

[pxx,f]=pwelch(S3,[],[],[],Fs);
plot(f,10*log10(pxx),'b','LineWidth',2)

ind = dsearchn(f,f0Vec);
scatter(f(ind),10*log10(pxxS(ind)),'r','filled');

xlim([0 80])
hold off


