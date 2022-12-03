% fractal analysis
% aim: separating fractal component >> calculate a in psd>>1/f^a
% currently work for single trial
% date '17-May-2022'
% calculate 'irasa' using fieldtrip toolbox


function [Fractalpowspctrm,Xpowspctrm,f] = AFfractalAnalysis(X,fs,L_Segment,Overlap,Freq_Lim)

t=0:1/fs:(length(X)/fs)-(1/fs);

addpath E:\CurrentWork\NMP_Wiener\220426\fieldtrip-20220514
% addpath H:\CurrentWork\NMP_Wiener\1\Function\fieldtrip-20220514

% addpath H:\NMP_Wiener\220426\fieldtrip-20220514
% addpath G:\NMP_Wiener\220426\fieldtrip-20220514
ft_defaults

% data
data.trial{1,1} = X';
data.time{1,1}  = t;
data.label{1}     = 'chan';
data.trialinfo(1,1) = 1;

% chunk into segments
cfg               = [];
cfg.feedback	='no';
cfg.verbose	='no';
cfg.length        = L_Segment;
cfg.overlap       = Overlap;
data              = ft_redefinetrial(cfg, data);

% compute the fractal and original spectra
cfg               = [];
cfg.foilim        = Freq_Lim;        %example: [2 200] 
cfg.pad           = 'nextpow2';
cfg.method        = 'irasa';
cfg.output        = 'fractal';

fractal = ft_freqanalysis(cfg, data);
cfg.output        = 'original';
original = ft_freqanalysis(cfg, data);

rmpath E:\CurrentWork\NMP_Wiener\220426\fieldtrip-20220514
% rmpath H:\CurrentWork\NMP_Wiener\1\Function\fieldtrip-20220514
% rmpath H:\NMP_Wiener\220426\fieldtrip-20220514
% rmpath G:\NMP_Wiener\220426\fieldtrip-20220514

f= fractal.freq;
Xpowspctrm =original.powspctrm;
Fractalpowspctrm =fractal.powspctrm;
end