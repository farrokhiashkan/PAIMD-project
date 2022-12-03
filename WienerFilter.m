function [ObsT_Filter,ObsT_Filter2,b] = WienerFilter(Obs,fwhm,peakf,SpecOption,f0Vec,Width,Fractal_fit)


%% F0:
fs=SpecOption.Fs;
fwhm = ceil(fwhm);

%% Plot Signal

[Pxx,F]=periodogram(Obs,[],[],fs,'power');
PXX2=log10(Pxx(2:end,:));
F = F(2:end);


%%%%%%%%%%%%
Fstart =dsearchn(F,peakf-fwhm);
Fend = dsearchn(F,peakf+fwhm);
PXX3=PXX2;
PXX3(Fstart:Fend)=[];
F2=F;
F2(Fstart:Fend)=[];
%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%
Ind = dsearchn(f0Vec',peakf);
f0Vec(Ind)=[];
Width(Ind)=[];
Width=Width*2;
for i=1:length(f0Vec)
%     f11 =dsearchn(F2,f0Vec(i)-3);
%     f22 = dsearchn(F2,f0Vec(i)+3);

    f11 =dsearchn(F2,f0Vec(i)-Width(i));
    f22 = dsearchn(F2,f0Vec(i)+Width(i));
    
    
    PXX3(f11:f22)=[];%mean(mean(PXX3(f11-1:f11))+mean(PXX3(f22:f22+1)));
    F2(f11:f22)=[];
end
%%%%%%%%%%%%%%

PXX3 = interp1(F2,PXX3,F,'pchip');

Fstart =dsearchn(F2,SpecOption.Freq_Lim(1));
Fend = dsearchn(F2,SpecOption.Freq_Lim(2));

Xpred = [ones(length(F2(Fstart:Fend)),1) log10(F2(Fstart:Fend))];
b(1,:) = lscov(Xpred,PXX3(Fstart:Fend,1));

% L_Segment = SpecOption.L_Segment;
% Overlap = SpecOption.Overlap;
% Freq_Lim = SpecOption.Freq_Lim;
% [Fractalpowspctrm,~,F] = AFfractalAnalysis(zscore(Obs),fs,L_Segment,Overlap,Freq_Lim);
%    
% Xpred = [ones(length(F),1) log10(F')];
% b(1,:) = lscov(Xpred,log10(Fractalpowspctrm'));
    

%% Design Filter:

N = length(Obs);

ObsF = fft(Obs,N) ;
PFObs = abs(ObsF).^2;

Freq = (-N/2:N/2-1)*(fs/N); % zero-centered frequency range
Freq = fftshift(Freq);
FreqL = Freq(1:N/2);


B = (log10(PFObs(1:N/2)))' - (b(1,2)*log10(FreqL(1:end,:)));
b(1,1) = mean(B(2:end));

Snn = b(1,1) + ((b(1,2))*log10(FreqL(1:end,:)));% + X;
st = b(1,1) + b(1,2)*log10((N/2)*(fs/N));% + X;
Snn2 = [st flip(Snn(2:end))];% + X;
Snn=fftshift([Snn Snn2]);
Snn = 10.^Snn;

if ~isnan(Fractal_fit)
    FracF_Phase = angle(fft(Fractal_fit,N)) ;
    Snn = Snn.*exp(1i*fftshift((FracF_Phase)));
end





%=====
Freq = fftshift(Freq);
PFObs = fftshift(PFObs);
F0=peakf;
F_center = dsearchn(Freq',0);
L = [F_center-round(F0*(N/fs)) F_center+round(F0*(N/fs))];

HBW = fwhm;
W = round(HBW/((fs)/N));

if ~isnan(Fractal_fit)
    R0(:,1)=PFObs.*exp(1i*angle(fftshift(ObsF)));
else
    R0(:,1)=PFObs;
end


R0(L(1)-W:L(1)+W,1)=Snn(L(1)-W:L(1)+W);
R0(L(2)-W:L(2)+W,1)=Snn(L(2)-W:L(2)+W);
R0 = fftshift(R0);
PFObs = fftshift(PFObs);
Snn=R0;

%=== Filter:
if isnan(Fractal_fit)
    Sss = PFObs-Snn;
else
    Sss = abs(PFObs.*exp(1i*angle(fftshift(ObsF))) - Snn);
end

% Sss  = PFObs.*exp(1i*phase(ObsF)) - Snn;

Wfilter = (sign(Sss).*abs(Sss))./(abs(Sss)+abs(Snn));
% Wfilter = (abs(PFObs-Snn)./(abs(PFObs)));
Wfilter2 = (sign(Snn).*abs(Snn))./(abs(Snn)+abs(Sss));

% Wfilter = (Sss./(Snn+Sss)) ;

%=== Denoising Signal
ObsF_Filter = ObsF.*Wfilter ;
ObsT_Filter = real(ifft(ObsF_Filter));

ObsF_Filter2 = ObsF.*Wfilter2 ;
ObsT_Filter2 = real(ifft(ObsF_Filter2));


end
