%% this code calculates the Peaks in spectrum of signals
%% you need field trip to use this code >> IRASA calculated as part of peak_finder algorithm
%% Ashkan Farrokhi
%% IUST
%% Date: '29-Nov-2022'
%% Describtion:

%= this vesion of code only work on averaged signal: at final step the code
% find peak in mean of whitened signal (across trials) 

%===Input:
%           Signal: [sample * trial(or channel)] 
%           Fs: frequency sampling 
%           NumberOfComp: number of frequncy components (F0) that returned

%===Output:
%           F0_estimate: average frequency of components


%% Function:
function [F0_estimate,Width_estimate,pointer]=AFsearchForF0(Signal,Fs,NumberOfComp,MaxPointer,SpecOption,Fvec)

L_Segment = SpecOption.L_Segment;
Overlap = SpecOption.Overlap;
Freq_Lim = SpecOption.Freq_Lim;


b = nan(size(Signal,2),2);
x2 = nan(size(Signal));
for i =1:size(Signal,2)
    
    % calculate 'irasa' using fieldtrip toolbox
    [Fractalpowspctrm,Xpowspctrm(:,i),F] = AFfractalAnalysis(Signal(:,i),Fs,L_Segment,Overlap,Freq_Lim);
    
    %Fractalpowspctrm = smoothdata(Fractalpowspctrm,'gaussian',100);    
    
    Xpred = [ones(length(F),1) log10(F')];
    b(i,:) = lscov(Xpred,log10(Fractalpowspctrm'));
    
    h=.001;
    a=((-b(i,2))/2);
    S = Signal(:,i);
    n  = numel(S);
    J  = 0:(n-1);
    G1 = gamma( J+1 );
    G2 = gamma( a+1-J );
    s  = (-1) .^ J;
    M  = tril( ones(n) );
    R  = toeplitz( S(:)' );
    T  = meshgrid( (gamma(a+1)/(h^a)) * s ./ (G1.*G2) );
    x2(:,i)  = reshape(sum( R .* M .* T, 2 ), size(S));
end
Xpowspctrm = mean(Xpowspctrm,2);

% figure
[sst,f]=wsst(mean(x2,2),Fs,'amor');


Fstart = dsearchn(f',SpecOption.Freq_Lim(1));
Fend = dsearchn(f',SpecOption.Freq_Lim(2));

sst(1:Fstart-1,:)=0;
sst(Fend+1:end,:)=0;

[pks,lsor,Width,~] = findpeaks(Xpowspctrm,F,'SortStr','descend','Annotate','extents','WidthReference','halfprom');

Comp=0;
pointer=NumberOfComp-1;
while Comp<MaxPointer
    pointer = pointer+1;
    [fridge,~] = wsstridge(sst,100,f,'NumRidges',pointer);
    
%    figure
%    plot(1:1:size(Signal,1),fridge,'k','linewidth',1); 
%    xrec = iwsst(sst,iridge);
    
    f0 = mean(fridge);
    
    FF =dsearchn(lsor',f0');
    Width0 = Width(FF);
    pks_current= pks(FF);
    FF= lsor(FF);

    
    [F0,ia,~] = unique(FF);
    pks_current= pks_current(ia);
    Width0 = Width0(ia);
    
    [~,IndPks]=sort(pks_current,'descend');
    Width0 = Width0(IndPks);
    F0 = F0(IndPks);

    SumCond=0;
    F0_estimate=nan(NumberOfComp,1);
    Width_estimate=nan(NumberOfComp,1);
    for i=1:length(F0)
        for Osc = 1:size(Fvec,1)
            Dist(1,1) = (sign(F0(i)-Fvec(Osc,1)));
            Dist(1,2) = (sign(F0(i)-Fvec(Osc,2)));
            Dist = sum(Dist,2);
            if Dist==0 
                SumCond= SumCond+1;
                F0_estimate(SumCond)=F0(i);
                Width_estimate(SumCond)=Width0(i);
            end
            
        end %Osc
    end
    
    if SumCond==NumberOfComp
        Comp=MaxPointer;
    end
    
    if pointer>=MaxPointer
        Comp=MaxPointer; 
    end
end



end




