%% code for generating simulation signal
%% Signal = IMFs + Fractal
%% Ashkan Farrokhi
%% Date '29-Nov-2022'
%% IUST
%% Initialization:
clear;
clc; 
close all;

%% Signal spec:
Fs=1000;                                                                   % frequency sampling(Hz) of the signal
L_Trial=4000;                                                              % signal's length (for each trial)
N_Trial=1;                                                                 % number of trials
t=0:1/Fs:(L_Trial/Fs)-(1/Fs);                                              % time points

Snr_WhiteNoise = Inf; 
%Note: this SNR is different from what defined in paper
%% Generate Brownian Noise: (fractal component spec)
% for noise with different exponent changed the Color parameter
cn = dsp.ColoredNoise('Color','brown','SamplesPerFrame',L_Trial,'NumChannels',N_Trial);
X= cn();
%% periodic components (IMFs' spec): 
s1 = nan(L_Trial,N_Trial);
s2 = nan(L_Trial,N_Trial);
s3 = nan(L_Trial,N_Trial);
R = nan(L_Trial,N_Trial);

%define the frequency of oscillations
f1=10; %Hz
f2=20; %Hz
f3=40; %Hz
%% Generate Periodic Signals
for i=1:N_Trial
    
    %=== s1
    a1=5;
    a=a1+(.1*(a1*(2*rand(1)-1)));
    f=f1 +(2*rand(1)-1);
    s1(:,i)= a*cos(2*f*pi*t);
        
    %=== s2    
    a1=6;
    a=a1+(.1*(a1*(2*rand(1)-1)));
    f=f2 +(2*rand(1)-1);
    df1=(pi/2)*(2*rand(1)-1);
    df2=(pi/2)*(2*rand(1)-1);
    df3=0;
    s2(:,i)= a*(1+.4*sin(2*pi*(3)*t+df1)).*sin(2*f*pi*t+.6*sin(2*(4)*pi*t+df3)+df2);

    
    %=== s3
    a1=5;
    a=a1+(.1*(a1*(2*rand(1)-1)));
    f=f3+(2*rand(1)-1);
    df1=(pi/2)*(2*rand(1)-1);
    df2=(pi/2)*(2*rand(1)-1);
    df3=0;
    s3(:,i)= a*(1+.2*sin(2*pi*(5)*t+df1)).*sin(2*f*pi*t+.1*sin(2*(6)*pi*t+df3)+df2);
      
    R(:,i) = s1(:,i) + s3(:,i)+s2(:,i);
end

%% observed simulated signal >> (periodic + aperiodic + white noise)
R = awgn(R,Snr_WhiteNoise,'measured');                                     %add white noise
Signal = X+R;                                                              %add fractal background


%% Save results
Report.Signal.BrownianNoise=X;                                             % Background
Report.Signal.S1=s1;                                                       % Source 1
Report.Signal.S2=s2;                                                       % Source 2
Report.Signal.S3=s3;                                                       % Source 3
Report.Signal.SourcePlusWhiteNoise=R;                                      % S1+S2+S3+(White_Noise)
Report.Signal.Observation=Signal;                                          % R+X

Report.Propety.SNR_whiteNoise_onPeriodicPart=Snr_WhiteNoise;
Report.Propety.Fs=Fs;

save('LFP_Simulated_Signal.mat','Report')


%% plot 
figure(1)

subplot(6,1,1), plot(t,s1), ylabel('S1')
subplot(6,1,2), plot(t,s2), ylabel('S2')
subplot(6,1,3), plot(t,s3), ylabel('S3')

subplot(6,1,4), plot(t,R), ylabel('R')
subplot(6,1,5), plot(t,X), ylabel('X')
subplot(6,1,6), plot(t,Signal), ylabel('Signal')


figure(2), hold on
title('PSD')
[pxx,f]=pwelch(Signal,[],[],[],Fs);
plot(f,10*log10(pxx),'k','LineWidth',3)

[pxx,f]=pwelch(s1,[],[],[],Fs);
plot(f,10*log10(pxx),'b','LineWidth',2)

[pxx,f]=pwelch(s2,[],[],[],Fs);
plot(f,10*log10(pxx),'b','LineWidth',2)

[pxx,f]=pwelch(s3,[],[],[],Fs);
plot(f,10*log10(pxx),'b','LineWidth',2)

[pxx,f]=pwelch(X,[],[],[],Fs);
plot(f,10*log10(pxx),'--r','LineWidth',2)



xlim([0 80])
hold off
