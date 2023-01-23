%
% Script to run the Pred_OrLoc experiment at OHBA.
% 
% JY (Sep, 2022)
%

clearvars; close all; clc;
addpath(genpath(pwd));

global devMode
devMode = true; %false; %


%% =========================================================== 
%% ===              Ask for user's inputs                  ===
%% =========================================================== 
isMEG    = false;
freqRIFT = 1440;
[SubjectInfo, EnvInfo] = initExpt( 'PredOrLo', devMode, isMEG, freqRIFT );


%% =========================================================== 
%% ===              Define the experiment                  ===
%% ===========================================================
global wPtr scr stim ctl proc

% === Open screen ===
[wPtr, wRect] = PsychImaging('OpenWindow', EnvInfo.screenID, [],...
                    EnvInfo.rectDisplay, [], [], [], 8);
Screen('BlendFunction', wPtr, 'GL_SRC_ALPHA', 'GL_ONE_MINUS_SRC_ALPHA');

% === The screen ===
scr              = [];
scr.rect         = wRect;
scr.idx          = EnvInfo.screenID;
scr.black        = BlackIndex(wPtr);
scr.white        = WhiteIndex(wPtr);
scr.gray         = (scr.white + scr.black)./2;
scr.ifi     = Screen('GetFlipInterval', wPtr);
scr.size    = EnvInfo.sizeDisplay;
scr.bgColor = repmat(scr.gray, [1,3]);
scr.viewDist= EnvInfo.viewDistance; %in cm
% double check IFI
if round(1 / scr.ifi) ~= 120
    warning('The IFI suggests a refresh rate of something else than 120 Hz!');
end
% the center pixel
[scr.xCtr, scr.yCtr] = RectCenter(scr.rect);

scr.yCtr = scr.yCtr - 0.2*scr.yCtr;

% the scale (pixels per degree)
cfg = []; cfg = scr; cfg.degrees = [1,1];
scr.ppdXY = JY_VisExptTools('deg2Pixel',cfg);
scr.ppdX  = mean(scr.ppdXY);
scr.ppdY  = mean(scr.ppdXY);
                
% === Define the stimulus (display duration and location) ===
stim = struct('Fix',[], 'dur',[], 'loc',[]);

% Fixation
stim.Fix.color   = repmat(scr.black, [1,3]);
stim.Fix.size    = round(0.2 * scr.ppdX);
stim.Fix.rectFix = CenterRectOnPoint([0 0 stim.Fix.size stim.Fix.size], scr.xCtr, scr.yCtr);            

% Dimension Cue
stim.dimcue.ovalrectpix    = ceil( [0 0 1.2 1.2] .* scr.ppdX );
stim.dimcue.ovalrects(:,1) = CenterRectOnPoint(stim.dimcue.ovalrectpix , scr.xCtr - 1.75*scr.ppdX, scr.yCtr);
stim.dimcue.ovalrects(:,2) = CenterRectOnPoint(stim.dimcue.ovalrectpix , scr.xCtr + 1.75*scr.ppdX, scr.yCtr);

% Spatial Cue
tmpCue4Valid = {'R','D'}; %rectangle vs. diamond
stim.spacue.cue4valid = tmpCue4Valid{ mod( str2double(SubjectInfo.SubjectID),2) +1 };
stim.spacue.rtg_rectL = stim.dimcue.ovalrects(:,1);
stim.spacue.rtg_rectR = stim.dimcue.ovalrects(:,2);
[xR,~]  = RectCenter(stim.spacue.rtg_rectR);
[xL,~]  = RectCenter(stim.spacue.rtg_rectL);
radius  = floor( stim.dimcue.ovalrectpix(3) ./ 2 .* sqrt(2) );
stim.spacue.dm_rectR(:,1) = [xR, xR-radius, xR, xR+radius];%x-coordinates
stim.spacue.dm_rectR(:,2) = [scr.yCtr+radius, scr.yCtr, scr.yCtr-radius, scr.yCtr]; %y-coordinates
stim.spacue.dm_rectL(:,1) = [xL, xL-radius, xL, xL+radius];%x-coordinates
stim.spacue.dm_rectL(:,2) = [scr.yCtr+radius, scr.yCtr, scr.yCtr-radius, scr.yCtr]; %y-coordinates

% Duration
stim.dur.precue   = 0.4;
stim.dur.dimcue   = 0.2;
stim.dur.cue2cue  = 0.8;
stim.dur.spacue   = 0.2; %0.25;
stim.dur.cue2stim = 0.8;
stim.dur.s1       = 0.1; %1/60 * 5;
stim.dur.s2       = 0.1; %1/60 * 5;
stim.dur.ISI      = 0.6; 
stim.dur.MaxResp  = 2.5;
stim.dur.ITI      = [0.8, 1.2]; % NEED CONFIRM

% Location and size
stim.loc.xCtrDeg  = 9.8; %from Zhigalov & Jensen (2020) Hum Brain Mapp
stim.loc.xCtrPix  = stim.loc.xCtrDeg * scr.ppdX;
stim.loc.yCtrDeg  = 9.8; %from Zhigalov & Jensen (2020) Hum Brain Mapp
stim.loc.yCtrPix  = stim.loc.yCtrDeg * scr.ppdY;
stim.loc.rectDeg  = [0 0 5.7 5.7]; %from Zhigalov & Jensen (2020) Hum Brain Mapp
stim.loc.rectPix  = ceil( stim.loc.rectDeg .* scr.ppdX );
stim.loc.rect_R   = CenterRectOnPoint(stim.loc.rectPix, scr.xCtr + stim.loc.xCtrPix, stim.loc.yCtrPix + scr.yCtr);
stim.loc.rect_L   = CenterRectOnPoint(stim.loc.rectPix, scr.xCtr - stim.loc.xCtrPix, stim.loc.yCtrPix + scr.yCtr);

% define my bandpass-filtered patch
stim.patch.sizedeg   = stim.loc.rectDeg(3);
stim.patch.sizepix   = ceil( stim.patch.sizedeg .* scr.ppdX );
stim.patch.freq_mean = 2.5; %in cycles/deg
stim.patch.freq_sd   = 0.5; %0.5; %stim.patch.freq_mean / 3; %JY: arbitrary
stim.patch.ori_kappa = 400;
stim.patch.ori_mean  = []; %titrated
stim.patch.patchlum  = 0.5000;
stim.patch.patchcon  = 0.5; %the contrast when RIFTing
stim.patch.gauss_mask= 0;
stim.patch.gauss_sd  = [];
stim.noise           = stim.patch;

% define my aperture 
m.maskSiz = stim.patch.sizepix;
m.outerR  = stim.patch.sizedeg./2;
m.ppd     = scr.ppdX;
m.degSmoo = 0.5; %in dva
mm        = JY_VisExptTools('make_smooth_circular_mask', m);
stim.aperture = 1-mm;

% === Define the response keys and etc. ===
ctl.keySame   = KbName('S');
ctl.keyDiff   = KbName('D');
ctl.DeviceNum = -1;
ctl.keyQuit   = KbName('Q');
ctl.keyValid  = [ctl.keySame, ctl.keyDiff, ctl.keyQuit];


% === Define the experimental procedure ===
proc                 = [];
proc.nTrialsPerBlock = 64;
proc.nBlocks         = 8;
proc.nTrialsTotal    = proc.nTrialsPerBlock * proc.nBlocks;

% Possible trial configurations: Task=Loc, Validity(Val/Neutral) x Ori(H/V) x Loc2(Same/Diff)
[Task, CueVal, S1Loc,S2Loc, S1Ori,S2Ori] = ndgrid( ['S','O'], ['V','N'], ['L','R'],['L','R'], ['H','V'],['H','V'] );
if mod( proc.nTrialsTotal, numel(S1Loc(:)) )~=0
    error( 'proc.nTrialsTotal has to be multiples of the number of possible trial config!');
else
    nRpts = proc.nTrialsTotal ./ numel(S1Loc(:));
end

% Shuffle the trials
idx = randperm( proc.nTrialsTotal );
tmp = {'Task','CueVal','S1Loc','S2Loc','S1Ori','S2Ori'};
for iV = 1:numel(tmp)
    eval( sprintf('proc.%s = repmat( %s(:), [nRpts,1] );', tmp{iV}, tmp{iV}) );
    eval( sprintf('proc.%s = proc.%s( idx );',tmp{iV}, tmp{iV}) );
end


% === Define the eyetracker ===
global ET
if EnvInfo.eyelink
    
    ET.use  = 1;
    ET.rad  = scr.ppdX * 2; %radius of allowable movement in pixels
    ET.cx   = scr.xCtr;
    ET.cy   = scr.yCtr;
    ET.dir  = 'eyedata';
    
    ET.file = sprintf('OrLo_%s', num2str(SubjectInfo.SubjectID));
    
    if ~exist(ET.dir, 'dir'); mkdir(ET.dir); end %make the directory if necessary
    
else
    
    ET.use = 0;
    
end


% ==== SAVE THE PARAMETERS =====
save( SubjectInfo.datafile );



%% =========================================================== 
%% ===               photodiode placement                  ===
%% ===========================================================
Screen('FillRect', wPtr, [repmat(scr.gray,[1,3]),1], scr.rect );
switch stim.spacue.cue4valid
    case 'R'
        DrawFormattedText(wPtr, 'Follow the rectangle!', 'center', 'center', scr.black);
    case 'D'
        DrawFormattedText(wPtr, 'Follow the diamond!', 'center', 'center', scr.black);
end
Screen('Flip', wPtr); pause(2); KbWait(ctl.DeviceNum); 



%% =========================================================== 
%% ===            calibrate the eye tracker                ===
%% ===========================================================
if ET.use
    % Initialize eye tracker
    [ET.el exitFlag] = rd_eyeLink('eyestart_customize', wPtr, ET.file);
    if exitFlag, return, end
    
    % Calibrate the eye tracker
    [cal exitFlag] = rd_eyeLink('calibrate_customize',wPtr, ET.el);
    if exitFlag, return, end
end



%% =========================================================== 
%% ===                 Run the experiment                  ===
%% =========================================================== 
results = cell( proc.nBlocks, proc.nTrialsPerBlock );
for b = 1:proc.nBlocks
    
    % === show block instruction ===
    Screen('FillRect', wPtr, [repmat(scr.gray,[1,3]),1], scr.rect );
    bStartText = ['Press any key to start Block ', num2str(b), '.'];
    DrawFormattedText(wPtr, bStartText, 'center', 'center', scr.black);
    Screen('Flip', wPtr); pause(2); KbWait(ctl.DeviceNum);
    
    % === start eye-tracking ===
    if ET.use, Eyelink('StartRecording'); end
    
    % === buffer ===
    Screen('FrameOval', wPtr, scr.black, stim.Fix.rectFix, 4, 4); %fixation dot
    bStart = Screen('Flip', wPtr);
    
    % === loop through the trials ===
    for ii = 1:proc.nTrialsPerBlock
        
        % Define the trial
        trialcfg = trialdef_PredOL_BEH( wPtr, b, ii );
        if ii==1
            trialcfg.tStart = bStart + 1 + trialcfg.ITI;
        else
            trialcfg.tStart = results{b,ii-1}.resp.tEnd + trialcfg.ITI;
        end
        
        % Save the trial configuration parameters
        results{b,ii}.cfg  = trialcfg;
        results{b,ii}.resp = [];
        
        % Show trial
        results{b,ii}.resp = trialfun_PredOL_BEH(wPtr, trialcfg);
        
        % End trial by removing everything but fixation
        Screen('FrameOval', wPtr, scr.black, stim.Fix.rectFix, 2, 2); %fixation dot
        Screen('Flip', wPtr);
    end
    
    % === stop eye-tracking ===
    if ET.use, Eyelink('StopRecording'); end
    
    % === end of the block instruction ===
    bEndText = ['End of Block ', num2str(b), '.\n\n',...
                'Take a break while we save your data...'];
    Screen('FillRect', wPtr, [repmat(scr.gray,[1,3]),1], scr.rect );  
    DrawFormattedText(wPtr, bEndText, 'center', 'center', scr.black);
    Screen('Flip', wPtr);
    
    % === save data ===
    if exist( SubjectInfo.datafile, 'file')
        eval( sprintf('results%s = {results{%s,:}};', num2str(b), num2str(b)) );
        eval( sprintf('save( SubjectInfo.datafile, ''results%s'', ''-append'' );', num2str(b)) );
    else
        save( SubjectInfo.datafile );
    end
    pause(3); 
    
end


%% =========================================================== 
%% ===                    Clean up                         ===
%% =========================================================== 
sca; 
if ET.use, rd_eyeLink('eyestop',wPtr,{ET.file, ET.dir}); end
if EnvInfo.sendtrig
    Datapixx('SetPropixxDlpSequenceProgram', 0); %Revert to standard 120Hz refresh rate
    Datapixx('RegWrRd');
    Datapixx('Close');
end
disp('All good! Do not forget to transfer data!');