%% PAIMD >> decompose signal into paeriodic and aperiodic components
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
X = Report.Signal.Observation;                                             % Sum(Sources)+Background
Fs=Report.Propety.Fs;
%% Spec 
SpecOption.Fs=1000;                                                        % sampling frequency
SpecOption.Freq_Lim=[2 400];                                               % frequency limitation for processing
%% Set center frequency>>> use for initial condition
% the center frequency could detected by peak_finder algorithm or set by user
f1= 10;     Width1= 1;   
f2= 20;     Width2= 1;
f3= 40;     Width3= 1;
f0Vec = [f1 f2 f3];
Width = [Width1 Width2 Width3];
%% S1: estimate the first Source:
[IMF_final(:,1),Bfractal_1,ErrorTotal2]=AFPAIMD(X,f0Vec(1),SpecOption,f0Vec,Width);
[Pearson(1),NRMSE(1),R2(1), ~, ~]= MH_GoodnessFit(S1,IMF_final(:,1),0,0);

clc; 
disp('Pearson')
disp(Pearson)
disp('NRMSE')
disp(NRMSE)
disp('R2')
disp(R2)

% plot: 
figure(1), sgtitle('S1')
subplot(2,1,1)
plot(S1(1:end),'k'), hold on
plot(IMF_final(1:end),'--b'), ylabel('signal')
hold off

subplot(2,1,2)
plot(angle(hilbert(S1(1:end))),'k') , hold on
plot(angle(hilbert(IMF_final(2:end))),'--b'), ylabel('Phase')
hold off
legend('Signal','Estimation')
drawnow 
%% estimate the second Source:
[IMF_final(:,2),Bfractal_2]=AFPAIMD(X,f0Vec(2),SpecOption,f0Vec,Width);
[Pearson(2),NRMSE(2),R2(2), ~, ~]= MH_GoodnessFit(S2,IMF_final(:,2),0,0);

clc; 
disp('Pearson')
disp(Pearson)
disp('NRMSE')
disp(NRMSE)
disp('R2')
disp(R2)

% plot: 
figure(2), sgtitle('S2')
subplot(2,1,1)
plot(S2(1:end),'k'), hold on
plot(IMF_final(:,2),'--b'), ylabel('signal')
hold off

subplot(2,1,2)
plot(angle(hilbert(S2(1:end))),'k') , hold on
plot(angle(hilbert(IMF_final(:,2))),'--b'), ylabel('Phase')
hold off
legend('Signal','Estimation')
drawnow 
%% estimate the third Source:
[IMF_final(:,3),Bfractal_3]=AFPAIMD(X,f0Vec(3),SpecOption,f0Vec,Width);
[Pearson(3),NRMSE(3),R2(3), ~, ~]= MH_GoodnessFit(S3,IMF_final(:,3),0,0);

clc; 
disp('Pearson')
disp(Pearson)
disp('NRMSE')
disp(NRMSE)
disp('R2')
disp(R2)

% plot: 
figure(3), sgtitle('S3')
subplot(2,1,1)
plot(S3(1:end),'k'), hold on
plot(IMF_final(:,3),'--b'), ylabel('signal')
hold off

subplot(2,1,2)
plot(angle(hilbert(S3(1:end))),'k') , hold on
plot(angle(hilbert(IMF_final(:,3))),'--b'), ylabel('Phase')
hold off
legend('Signal','Estimation')
drawnow 
%% Background:
Fractal = X - sum(IMF_final,2);
True_Bractal = Report.Signal.BrownianNoise;


%% Plot PSD
figure(4)
hold on

[PX,F]=pwelch(X,[],[],[],Fs);figure(gcf)
plot(F,10*log10(PX),'k','LineWidth',2)

[PTF,F]=pwelch(True_Bractal,[],[],[],Fs);figure(gcf)
plot(F,10*log10(PTF),'r','LineWidth',2)


[PImf,F]=pwelch(IMF_final,[],[],[],Fs);figure(gcf)
plot(F,10*log10(PImf),'b','LineWidth',1.2)

[Pfractal,F]=pwelch(Fractal,[],[],[],Fs);figure(gcf)
plot(F,10*log10(Pfractal),'--g','LineWidth',1.2)

xlim([0 70])
ylim([-25 inf])
hold off
legend('Signal','Fractal','Estimated IMFs','','', 'Estimated Fractal activity')