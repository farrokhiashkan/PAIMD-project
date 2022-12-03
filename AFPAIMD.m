%% PAIMD algorithm >> decompose signal into paeriodic and aperiodic components
%% Ashkan Farrokhi
%% Date: '29-Nov-2022'
%% IUST
%% Input/Output
%=Input:
% X: input signal (sample*1)
% f0: initial condition for center frequency
% SpecOption: contain FS and frequency limit for processing>> for more information see Example 
% f0Vec: vector includes center frequency estimations
% width: estimation of bandwidth for estimated frequency

%=Output:
% IMF_final: estimated IMF
% fractal component (exponent and offset) estimation

%% Function:
function [IMF_final,Bfractal,ErrorTotal2]=AFPAIMD(X,f0,SpecOption,f0Vec,Width)
PlotStat=0;
G=0;
fs=SpecOption.Fs;
N=length(X);

% adding to signal to have even length for faster fft
if rem(N,2)==1
    X=[X;X(end)];
    N=N+1;
end %rem

t=0:1/fs:(N-1)/fs;
theta_Init= 2*pi*(f0+1)*t;


AlphaGrid =[.05, .1, .15, .2, .25, .3, .35, .4, .45 .5];

IMF_Total = nan(length(AlphaGrid),length(AlphaGrid),length(X));
ErrorTotal= nan(length(AlphaGrid),length(AlphaGrid));
ErrorTotal2= nan(length(AlphaGrid),length(AlphaGrid));
for alpha1=1:length(AlphaGrid)
    %alpha1
    for alpha2=1:length(AlphaGrid)
        theta = theta_Init;
        %alpha2
        Betarange=[0,.1,.2,.5, 1];
        Error=inf;
        Error2=inf;
        Cond=1;
        
        PointerError=0;
        pointer=0;
        IMF_Memory=[];
        Err2=[];
        Err1=[];
        while Cond==1
            pointer=pointer+1;
            
            %Interpolate the data to the theta space
            theta_fit = linspace(theta(1),theta(end),N);
            
            L = ceil((theta_fit(end)-theta_fit(1)) / (2*1*pi));
            L = L*fs/length(X);
            
            X_fit = interp1(theta,X,theta_fit,'pchip');
            
            % FFT
            Freq = (-N/2:N/2-1)*(fs/N); % zero-centered frequency range
            
            %%%%%%%%%%%% wiener filter
            Fractal_fit=nan;
            %  Osci_fit=nan;
            [X_fit,X_fit_Noise,Bfractal(alpha1,alpha2,:)] = WienerFilter(X_fit',AlphaGrid(alpha1)*L,L,SpecOption,f0Vec,Width,Fractal_fit);
            X_fit = X_fit';
            X_fit_Noise = X_fit_Noise';
            
            X_fit2 = (interp1(theta_fit,X_fit_Noise,theta,'pchip'));
            %  X_fit2=0;
            
            fft_X = fft(X_fit);
            fft_X = fftshift(fft_X);
            
            F_center = dsearchn(Freq',0);
            % L = ceil((theta_fit(end)-theta_fit(1)) / (2*1*pi));
            L = dsearchn(Freq',L);
            L = L-F_center;
            
            Shift_Right = circshift(fft_X,L);
            Shift_Left = circshift(fft_X,-L);
            
            LP_filter=.5*(1+ cos((pi*Freq)/(L)));
            LP_filter = ones(1,length(LP_filter));
            
            Lambda=AlphaGrid(alpha1);
            Ind_Edge = floor(Lambda*L)+0;
            Ind_Edge = [F_center-Ind_Edge, F_center+Ind_Edge];
            LP_filter([1:Ind_Edge(1),Ind_Edge(2):length(Freq)])=0;
            
            
            a = 1*(Shift_Left+Shift_Right).*LP_filter   ;
            a = real(ifft(ifftshift(a)));
            
            
            b = 1*(1i*(Shift_Left-Shift_Right)).*LP_filter   ;
            b = real(ifft(ifftshift(b)));
            
            
            at = 1*(interp1(theta_fit,a,theta,'pchip'));
            bt = 1*(interp1(theta_fit,b,theta,'pchip'));
            IMF = 1*(((at.*cos(theta))+(bt.*sin(theta))));
            
            
            %=== update Theta
            d_at = diff([at(1), at]);    d_at(1) = d_at(2);
            d_bt = diff([bt(1), bt]);    d_bt(1) = d_bt(2);
            
            dtheta = ((at.*d_bt) - (bt.*d_at)) ./ (at.^2 + bt.^2);
            dtheta(isnan(dtheta))=0;
            dtheta_fit = interp1(theta,dtheta,theta_fit,'pchip');
            
            fft_dtheta_fit = fft(dtheta_fit);
            fft_dtheta_fit = fftshift(fft_dtheta_fit);
            
            LP_filter=.5*(1+ cos((pi*Freq)/(L)));
            LP_filter = ones(1,length(LP_filter));
            
            Lambda=AlphaGrid(alpha2);
            Ind_Edge = floor(Lambda*L)+0;
            Ind_Edge = [F_center-Ind_Edge, F_center+Ind_Edge];
            LP_filter([1:Ind_Edge(1),Ind_Edge(2):length(Freq)])=0;
            
            dtheta = real(ifft(ifftshift(fft_dtheta_fit.*LP_filter)));
            theta_add = cumtrapz(t,dtheta);
            
            Beta=Betarange;%.1:1:100;
            theta_new = theta- Beta'.*theta_add ;
            Diff = diff(theta_new,[],2);
            
            Diff = sum(Diff<-10^-3,2);
            Ind = (Diff<=0);
            Ind  =find(Ind==1);
            Beta=Beta(Ind(end));
            theta_new = theta- Beta*theta_add ;
            
            
            if PlotStat==1
                figure(2)
                subplot(2,1,1)
                plot(t(G+1:end-G),Sig(G+1:end-G),'b'), hold on
%                                                 plot(t(G+1:end-G),s2,'r')
%                                                 plot(t(G+1:end-G),s3,'g')
                plot(t(G+1:end-G),IMF(G+1:end-G),'--k')
                hold off
                
                subplot(2,1,2)
                plot(phase(hilbert(Sig(G+1:end-G))),'b') , hold on
                %                                 plot(phase(hilbert(s2)),'r')
                %                                 plot(phase(hilbert(s3)),'g')
                plot(theta_new,'k')
                hold off
            end
            
            X2=X;%Sig;%X;%X-X_fit_Noise';
            [~,ErrorNew,~, ~, ~]= MH_GoodnessFit(X2(G+1:end-G),IMF(G+1:end-G)',0,0);
            
            X2=X-X_fit2';%%-X_fit_Noise';
            [~,ErrorNew2,~, ~, ~]= MH_GoodnessFit(X2(G+1:end-G),IMF(G+1:end-G)',0,0);
            
            if ErrorNew2>Error2
                PointerError=PointerError+1;
            else
                PointerError=0;
            end
            
            if pointer==400 || Beta==0 || norm(theta_add)<10^-10 || PointerError==10%|| norm(Beta*theta_add)/norm(theta)<.0001
                Cond=0;
                if pointer==1
                    IMFEnd = IMF;
                end
            else
                IMFEnd = IMF;
            end
            theta = theta_new;
            %                 thetaFfractal = phase(fft(X_fit2));  %%Ntry
            Error2=ErrorNew2;
            
            
            if pointer<=10
                IMF_Memory(pointer,:) = IMF;
                Err2(pointer) = ErrorNew2;
                Err1(pointer) = ErrorNew;
            else
                IMF_Memory(1:9,:) = IMF_Memory(2:10,:);
                IMF_Memory(10,:) = IMF;
                %
                Err2(1:9) = Err2(2:10);
                Err2(10) = ErrorNew2;
                %
                Err1(1:9) = Err1(2:10);
                Err1(10) = ErrorNew;
            end
            
            if PointerError==10
                IMFEnd = IMF_Memory(1,:);
                ErrorNew2 = Err2(1);
                ErrorNew = Err1(1);
            end
            Error=ErrorNew;
            
        end %Cond
        
        ErrorTotal(alpha1,alpha2) = Error;
        ErrorTotal2(alpha1,alpha2)=ErrorNew2;
        IMF_Total(alpha1,alpha2,:)= IMFEnd;
% %                     clc
% %                     disp('Ideal Correlation')
% %                     disp('Ideal NRMSE');
% %                     disp('Ideal R2')
% %                     disp(ErrorTotal);
% %                     disp('Correlation')
% %                     disp('NRMSE');
% %                     disp('R2')
% %                     disp(ErrorTotal2);
        %         catch
        %do nothing
    end %alpha2
end %alpha1

%%
I=min(ErrorTotal2,[],'all');
% min(ErrorTotal2,[],'all')
[I1,I2] = find(ErrorTotal2==I);

IMF_final = squeeze(IMF_Total(I1(1),I2(1),:));
Bfractal = Bfractal(I1(1),I2(1),:);
% disp([AlphaGrid(I1(1)) AlphaGrid(I2(1))])           
% ErrorTotal(I1,I2)

% close all
% 
% figure(20)
% subplot(2,1,1)
% plot(t(G+1:end-G),Sig(G+1:end-G),'b'), hold on
% plot(t(G+1:end-G),IMF_final(G+1:end-G),'--k')
% hold off
% 
% subplot(2,1,2)
% plot(phase(hilbert(Sig(G+1:end-G))),'b') , hold on
% plot(phase(hilbert(IMF_final)),'k')
% hold off
% 
% 
% figure
% [PX,F]=pwelch(X,[],[],[],256);figure(gcf)
% [PImf,F]=pwelch(IMF_final,[],[],[],256);figure(gcf)
% % [PSig,F]=pwelch(Sig,[],[],[],1000);figure(gcf)
% [PXImf,F]=pwelch(X-IMF_final,[],[],[],256);figure(gcf)
% % [PNoise,F]=pwelch(Noise,[],[],[],1000);figure(gcf)
% plot(F,10*log10(PX),'k')
% hold on
% plot(F,10*log10(PImf),'b')
% plot(F,10*log10(PXImf),'--r')
% % plot(F,10*log10(PNoise),'g')
% % plot(F,10*log10(PSig),'c')
% hold off


%
% disp('%%%%%%%%%%%%%%%%%%%%%%%%')
% [Pearson,NRMSE,R2, ~, ~]= my_GoodnessFit(X,IMF_final,0,0);


end