%
% Compare the photodiode's recorded signal with the RIFTing signal used by
% PsychToolbox.
% 
% And the ERFs to the task stimuli.
% 
% JY (Dec 2022)

clearvars; close all; clc;
addpath( '/home/zhouj/Documents/fieldtrip-20220819' );
ft_defaults;

set(0, 'defaultFigureRenderer', 'painter');

%% Behavioral file
beh = dir( ['/ohba/pi/knobre/joeyzhou/Pred_OrLo/Pilot/sub101/*.mat'] );
load( fullfile( beh.folder, beh.name), 'triggers','stim' );


%% Epoch the .fif file
fif = dir( ['/ohba/pi/knobre/joeyzhou/Pred_OrLo/Pilot/sub101/preproc/*/*preproc_raw.fif'] );
hdr = ft_read_header( fullfile( fif(1).folder, fif(1).name ) );
dat = [];
for iF = 1:numel(fif)
    disp( fif(iF).name );
    tmpdat = ft_read_data( fullfile( fif(iF).folder, fif(iF).name ) );
    dat    = [dat, tmpdat];
end

% ====== Define the channels of interest ======
iMEG    = cellfun(@(x) contains(x,'MEG'), hdr.label, 'uniformoutput',true);
iTrigg  = cellfun(@(x) contains(x,'STI101'), hdr.label, 'uniformoutput',true);

% ====== The Nth sample where we have the trigger ===
idxISI  = strfind( dat(iTrigg,:), [0 triggers.ISI1 triggers.ISI1] ) + 1; %JY: hard-coded 8
nTrials = numel(idxISI);

% ========== Segment (in samples) ===========
pre_dur      = 2;
post_dur     = 4;
pre_samples  = hdr.Fs * pre_dur;
post_samples = hdr.Fs * post_dur; 

% ========= Create fieldtrip data structure ==========
data0           = [];
data0.time      = cell( nTrials,1);
data0.trial     = cell( nTrials,1);
data0.label     = hdr.label;
data0.dimord    = 'chan_time';
data0.grad      = hdr.grad;
for iT = 1:nTrials
    sel = [-pre_samples:post_samples] + idxISI(iT);
    data0.sampleinfo(iT,:) = [-pre_samples, post_samples];
    data0.time{iT}  = -pre_dur:(1/hdr.Fs):post_dur;
    data0.trial{iT} = dat( :, sel );
end

%%
% ============= fourth-order Butterworth zero-phase filter ===========
cfg             = [];
cfg.bpfilter    = 'yes';
cfg.bpfreq      = [55, 75];
cfg.bpfilttype  = 'firws';
data            = ft_preprocessing( cfg, data0 );
cfg.hilbert     = 'angle';
datPhase        = ft_preprocessing( cfg, data0 );


%% ================ Visualize the trials ==============
%{
idxMEG = find( iMEG );
for iChan = 1:sum(iMEG)
% plotFFT( data0.time{1}, data0.trial{1}(idxMEG(iChan),:), '-'); hold on,
plotFFT( data.time{1}, data.trial{1}(idxMEG(iChan),:), '-'); hold on,
end

%%
for iT = 1:nTrials
    plotFFT( data0.time{iT}, data0.trial{iT}(iDiodeA,:), '-'); hold on,
end
%}


%% ================= Visualize the "evoked activity" ===============
%{
cfg              = [];
cfg.keeptrials   = 'yes';
D0               = ft_timelockanalysis(cfg, data0);

cfg              = [];
cfg.baseline     = [-0.5, 0];
cfg.baselinetype = 'rel';
ERPs             = ft_timelockbaseline( cfg, D0 );%

ERPs = ft_timelockanalysis( [], ERPs );

cfg              = [];
cfg.parameter    = 'avg';
cfg.xlim         = [0:0.05:0.2];
cfg.zlim         = [-1,1] * 3 * 1e-12;%'zeromax';
cfg.layout       = 'neuromag306planar_helmet.mat';
cfg.colorbar     = 'yes';
cfg.comment      = 'xlim';
cfg.commentpos   = 'title';
ft_topoplotER( cfg, ERPs );
set(gcf,'units','centimeters','position',[0 0 24 24]);
% print('RIFT_ERFs','-dpng','-r800');
%}


%% ======== time-freq analysis ==========
%{
cfg           = [];
cfg.output    = 'pow';
cfg.channel   = {'MEG'};
cfg.method    = 'mtmconvol';
cfg.taper     = 'hanning';
cfg.foi       = 2:1:30;
cfg.t_ftimwin = 3./cfg.foi;
cfg.toi       = -0.5:0.05:1.5;
TFR           = ft_freqanalysis(cfg, data0);
%}



%% Compute PLV (lags between -200 to 200 ms, in 1 ms steps)
iDiodeA = cellfun(@(x) contains(x,'MISC011'), datPhase.label, 'uniformoutput',true);
iDiodeB = cellfun(@(x) contains(x,'MISC012'), datPhase.label, 'uniformoutput',true);

% Compute instantaneous phase
lags_sec    = [-0.1:(1/hdr.Fs):0.2]; %[-0.2:0.001:0.2];
lags_sample = lags_sec .* hdr.Fs;

idxMEG      = find(iMEG);
timwin      = [0, stim.dur.ISI+stim.dur.s1];
[~,iBeg]    = min( abs(timwin(1)-data.time{1}) );
[~,iEnd]    = min( abs(timwin(end)-data.time{1}) );
selDiode    = iBeg:iEnd;
selMEG0     = (iBeg+lags_sample(1)):(iEnd+lags_sample(end));

inst_phase_diodeA = nan(nTrials, 1, numel(selDiode));
inst_phase_diodeB = nan(nTrials, 1, numel(selDiode));
inst_phase_MEG    = nan(nTrials, sum(iMEG), numel(selMEG0));
for iT = 1:nTrials
    
    % ========= the diode signal =========
    %{
    yDiodeA = data.trial{iT}(iDiodeA, selDiode);
    yDiodeB = data.trial{iT}(iDiodeB, selDiode);
     
    zDiode = hilbert( yDiodeA );
    inst_phase_diodeA(iT,:) = angle(zDiode);
    
    zDiode = hilbert( yDiodeB );
    inst_phase_diodeB(iT,:) = angle(zDiode);
    %}
    inst_phase_diodeA(iT,1,:) = datPhase.trial{iT}(iDiodeA,selDiode);
    inst_phase_diodeB(iT,1,:) = datPhase.trial{iT}(iDiodeB,selDiode);
    
    % ========= the MEG signal =========
    %{
    for iC = 1:sum(iMEG)
        datMEG = data.trial{iT}(idxMEG(iC), selMEG0);
        datMEG = hilbert( datMEG );
        inst_phase_MEG(iT,iC,:) = angle( datMEG );
    end
    %}
    inst_phase_MEG(iT,:,:) = datPhase.trial{iT}(iMEG,selMEG0);
end
    
%% compute PLV between the MEG and Diode =======
PLV0_DiodeA = nan( sum(iMEG), numel(lags_sample));
PLV0_DiodeB = nan( sum(iMEG), numel(lags_sample));
n_samples   = numel( selDiode );
for iL = 1:numel(lags_sample)
    selMEG = [ 1:numel(selDiode) ] + iL - 1;
    disp( [selMEG(1), selMEG(end)] );
    for iC = 1:sum(iMEG)
        PLV0_DiodeA(iC,iL) = mean(  abs( sum(exp(1i*(inst_phase_diodeA-inst_phase_MEG(:,iC,selMEG))),3) ) / n_samples   );
        PLV0_DiodeB(iC,iL) = mean(  abs( sum(exp(1i*(inst_phase_diodeB-inst_phase_MEG(:,iC,selMEG))),3) ) / n_samples   );
    end
end


%% Visualize the PLV
PLV_struct        = [];
PLV_struct.time   = lags_sec;
PLV_struct.dimord = 'chan_time';
PLV_struct.grad   = data.grad;
PLV_struct.label  = hdr.label( iMEG );
PLVa = PLV_struct;
PLVa.avg = PLV0_DiodeA; %squeeze( mean(PLV0_DiodeA,1) );
PLVb = PLV_struct;
PLVb.avg = PLV0_DiodeB; %squeeze( mean(PLV0_DiodeB,1) );

cfg                  = [];
cfg.parameter        = 'avg';
cfg.xlim             = [-0.1:0.05:0.2];%[0:0.05:0.2];
cfg.zlim             = [0.1, 0.25]; %
cfg.layout           = 'neuromag306planar_helmet.mat';
cfg.highlight        = 'on';
cfg.highlightcolor   = 'r';
cfg.highlightsymbol  = '.';
cfg.highlightsize    = 12;
cfg.baseline         = 'no';
cfg.colorbar         = 'yes';
cfg.comment          = 'xlim';
cfg.commentpos       = 'title';

save('PLV0.mat','PLVa','PLVb');

[~,chanIdxA]         = maxk( nanmean( PLVa.avg,2), 5);
cfg.highlightchannel = { PLVa.label(chanIdxA) };
ft_topoplotER( cfg, PLVa );
set(gcf,'units','centimeters','position',[0 0 36 24]);
print('RIFT_topoPLV_planar_DiodeRight','-dpng','-r800');

cfg.highlightcolor   = 'b';
[~,chanIdxB]         = maxk( nanmean( PLVb.avg,2), 5);
cfg.highlightchannel = { PLVa.label(chanIdxB) };
ft_topoplotER( cfg, PLVb );
set(gcf,'units','centimeters','position',[0 0 36 24]);
print('RIFT_topoPLV_planar_DiodeLeft','-dpng','-r800');

% cmbPLVa = ft_combineplanar([], PLVa);





%% =========== how do the ERFs look like in these RIFTy channels ========
% preISI = stim.dur.precue + stim.dur.dimcue + stim.dur.ISI + stim.dur.spacue;

cfg              = [];
cfg.keeptrials   = 'yes';
cfg.channel      = {'MEG'};
D0               = ft_timelockanalysis(cfg, data0);

cfg              = [];
cfg.baseline     = [-0.5, 0] - preISI;
cfg.baselinetype = 'rel';
ERPs             = ft_timelockbaseline( cfg, D0 );

% ERPs = ft_combineplanar( [], ERPs );
ERPs = ft_timelockanalysis( [], ERPs );

%%
close all; 
figure, clf, hold on,
set(gcf,'units','centimeters','position',[0 0 8 4]);
cfgfig         = [];
cfgfig.color   = {'r','b'};
cfgfig.chanIdx = {chanIdxA, chanIdxB};
cfgfig.stimdur = stim.dur;
cfgfig.ylim    = [-3, 4].*1e-12;
plotERFs( cfgfig, ERPs );
print('ERFs_rifty_channels','-dpng','-r800');

%%
figure, clf, hold on,
set(gcf,'units','centimeters','position',[0 0 8 4]);
preISI = stim.dur.precue + stim.dur.dimcue + stim.dur.ISI + stim.dur.spacue;
cfgfig.xlim    = [preISI-0.05, preISI + stim.dur.ISI1 + stim.dur.s1 + 0.3];
cfgfig.xticks  = cfgfig.xlim(1):0.1:cfgfig.xlim(2);
plotERFs( cfgfig, ERPs );
print('ERFs_rifty_channels_riftywin','-dpng','-r800');


%% ========= sub-function ============
function plotFFT( t, x, lnspec )

Fs= 1./diff(t(1:2));
t = t - Fs;
L = numel(t);


Y = fft(x);

P2 = abs(Y/L);
P1 = P2(1:L/2+1);
P1(2:end-1) = 2*P1(2:end-1);

f = Fs*(0:(L/2))/L;

f2plot = [0, 100];
idxBeg = 1;
idxEnd = find( f<= f2plot(2), 1, 'last' );

plot(f(idxBeg:idxEnd),P1(idxBeg:idxEnd), lnspec) 
% title("Single-Sided Amplitude Spectrum of X(t)")
% xlabel("f (Hz)")
% ylabel("|P1(f)|")

set(gca,'YScale','log','XScale','log');
xlim( [10^0, 10^2] );

end

%% =========== sub-function ============
function plotERFs( cfgfig, ERPs )

preISI = cfgfig.stimdur.precue + cfgfig.stimdur.dimcue + cfgfig.stimdur.ISI + cfgfig.stimdur.spacue;

for iLn = 1:numel(cfgfig.chanIdx)
    
    y = nanmean( ERPs.avg( cfgfig.chanIdx{iLn}, :), 1);
    
    plot( ERPs.time+preISI, y, cfgfig.color{iLn} );
end

set(gca,'TickDir','out','fontname','arial','fontsize',6);
xlabel( 'Time (s)','fontweight','bold' );
ylabel( 'Event-Related Response','fontweight','bold' );

dimcueON = cfgfig.stimdur.precue;
spacueON = dimcueON + cfgfig.stimdur.dimcue + cfgfig.stimdur.ISI;
riftON   = spacueON + cfgfig.stimdur.spacue;
stim1ON  = riftON + cfgfig.stimdur.ISI1;
stim2ON  = stim1ON + cfgfig.stimdur.s1 + cfgfig.stimdur.ISI2;

stim1OFF = stim1ON + cfgfig.stimdur.s1;

% box marking the rifting period
riftybox = area( [riftON, stim1OFF], [cfgfig.ylim(2),cfgfig.ylim(2)], cfgfig.ylim(1) );
riftybox.FaceColor = [0.8,0.8,0.8];
riftybox.EdgeColor = 'none';
riftybox.FaceAlpha = 0.2;

% lines marking the onsets of events
xline( dimcueON, 'k:' );
xline( spacueON, 'k:' );
xline( riftON, 'k:');
xline( stim1ON, 'k:');

ylim( cfgfig.ylim );

cfgfig.xlim   = ft_getopt(cfgfig,'xlim',[-0.5, max(stim2ON)+0.5]);
cfgfig.xticks = ft_getopt(cfgfig,'xticks',cfgfig.xlim(1):0.5:cfgfig.xlim(2));
xlim( cfgfig.xlim );
xticks( cfgfig.xticks );


end