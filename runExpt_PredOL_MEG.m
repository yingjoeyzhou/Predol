%
% Script to run the Pred_LoCo experiment at OHBA.
% 
% JY (Sep, 2022)
%

clearvars; close all; clc;
addpath(genpath(pwd));

global devMode
devMode = false; %true; %


%% =========================================================== 
%% ===              Ask for user's inputs                  ===
%% =========================================================== 
isMEG    = true;
freqRIFT = 1440;
[SubjectInfo, EnvInfo] = initExpt( 'PredOrLo', devMode, isMEG, freqRIFT );

EnvInfo.eyelink = 0;


%% =========================================================== 
%% ===              Define the experiment                  ===
%% ===========================================================
global wPtr scr stim ctl proc


%% ============ non-MEG stuff ============
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
stim.spacue.rtg_rectR = stim.dimcue.ovalrects(:,1);
stim.spacue.rtg_rectL = stim.dimcue.ovalrects(:,2);
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
stim.dur.ITI      = [0.8, 1.2];

% Location and size
stim.loc.xCtrDeg  = 7.0; %3.8; %9.8; %from Zhigalov & Jensen (2020) Hum Brain Mapp
stim.loc.xCtrPix  = stim.loc.xCtrDeg * scr.ppdX;
stim.loc.yCtrDeg  = 7.0; %3.8; %9.8; %from Zhigalov & Jensen (2020) Hum Brain Mapp
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
stim.patch.patchcon  = 1; %full contrast
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
if devMode==false
    ctl.keySame   = KbName('1!'); %Blue
    ctl.keyDiff   = KbName('4$'); %Red
else
    ctl.keySame   = KbName('S'); %Blue
    ctl.keyDiff   = KbName('D'); %Red
end
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


%% ======= MEG-specific stuff ======
% === Define the rapid invisible frequency tagging (RIFT) ===
global RIFT

PredolRIFT = load( 'Predol_RIFT.mat', 's1','s2', 'dur');

RIFT          = [];
RIFT.totaldur = 1.6; %JY: hard-coded

if RIFT.totaldur==PredolRIFT.dur %check duration of the tagging signal
    nPair = size(PredolRIFT.s1, 2);
    nSel  = proc.nTrialsPerBlock/4; %JY: hard-coded
    sel   = randperm( nPair, 16 );
    
    RIFT.s1 = PredolRIFT.s1(:, sel); %a subset of the pre-generated tagging signals is selected 
    RIFT.s2 = PredolRIFT.s2(:, sel);
    
    [sigLeft, sigIdx] = ndgrid( [1,2], 1:nSel );
    nRpts   = proc.nTrialsTotal ./ numel(sigLeft);
    idx     = randperm( proc.nTrialsTotal );
    sigLeft = repmat( sigLeft, [1,nRpts] );
    sigIdx  = repmat( sigIdx, [1,nRpts] );
    RIFT.whichsig = sigLeft( idx ); %whether signal from s1 or s2 is put on the left?
    RIFT.indexsig = sigIdx( idx ); %which siganl from s1/s2 is used?
end

RIFT.nFrame0   = stim.dur.cue2stim .* 120; %JY: hard-coded 120 Hz, which is the graphic card's refresh rate
RIFT.nFrameS1  = stim.dur.s1 .* 120; %JY: hard-coded 120 Hz
RIFT.nFrameISI = stim.dur.ISI .* 120; %JY: hard-coded 120 Hz
RIFT.nFrameS2  = stim.dur.s2 .* 120; %JY: hard-coded 120 Hz
RIFT.nFrameTot = RIFT.nFrame0 + RIFT.nFrameS1 + RIFT.nFrameISI + RIFT.nFrameS2;
RIFT.frameOnset= 1 + [0, RIFT.nFrame0, RIFT.nFrame0+RIFT.nFrameS1,...
                        RIFT.nFrame0+RIFT.nFrameS1+RIFT.nFrameISI];
                        
RIFT.Fix      = stim.Fix;
RIFT.Fix.size = ceil( RIFT.Fix.size ./ 2 );

RIFT.rectPix  = stim.loc.rectPix ./ 2;
RIFT.xCtrQuad = [ scr.xCtr.*0.5, scr.xCtr.*1.5, scr.xCtr.*0.5, scr.xCtr.*1.5 ];
RIFT.yCtrQuad = [ scr.yCtr.*0.5, scr.yCtr.*0.5, scr.yCtr.*1.5, scr.yCtr.*1.5 ] - scr.yCtr.*0.2; %JY: hard-coded to move fixation cross up!

for iQuad = 1:4 %becuase the the full screen has been divided into 4
    RIFT.rect_L{iQuad} = CenterRectOnPoint(RIFT.rectPix, RIFT.xCtrQuad(iQuad)-stim.loc.xCtrPix./2, RIFT.yCtrQuad(iQuad)+stim.loc.yCtrPix./2);
    RIFT.rect_R{iQuad} = CenterRectOnPoint(RIFT.rectPix, RIFT.xCtrQuad(iQuad)+stim.loc.xCtrPix./2, RIFT.yCtrQuad(iQuad)+stim.loc.yCtrPix./2);
end
% ========= Define the patch for photodiode measurement =========
RIFT.patch           = [];
RIFT.patch.sizedeg   = 1;
RIFT.patch.sizepix   = ceil( RIFT.patch.sizedeg .* scr.ppdX ./ 2 );
RIFT.patch.rectPix   = [0 0 RIFT.patch.sizepix, RIFT.patch.sizepix];

% diode location on the right
RIFT.patch.xCtr      = [scr.xCtr, scr.rect(3), scr.xCtr, scr.rect(3)] - RIFT.patch.sizepix/2; %put it to the right bottom corner
RIFT.patch.yCtr      = [scr.yCtr, scr.yCtr, scr.rect(4), scr.rect(4)] - RIFT.patch.sizepix/2; %put it to the right bottom corner
RIFT.patch.locRectsR = cell(4, 1); %photodiode positions on the right
for q = 1:4
    RIFT.patch.locRectsR{q} = CenterRectOnPoint(RIFT.patch.rectPix, ...
                                    RIFT.patch.xCtr(q), RIFT.patch.yCtr(q));
end

% diode location on the left
RIFT.patch.xCtr      = [scr.xCtr, scr.rect(1), scr.xCtr, scr.rect(1)] + RIFT.patch.sizepix/2; %put it to the right bottom corner
RIFT.patch.yCtr      = [scr.yCtr, scr.yCtr, scr.rect(4), scr.rect(4)] - RIFT.patch.sizepix/2; %put it to the right bottom corner
RIFT.patch.locRectsL = cell(4, 1); %photodiode positions on the right
for q = 1:4
    RIFT.patch.locRectsL{q} = CenterRectOnPoint(RIFT.patch.rectPix, ...
                                    RIFT.patch.xCtr(q), RIFT.patch.yCtr(q));
end

% diode stimuli
mR.maskSiz = RIFT.patch.sizepix;
mR.outerR  = RIFT.patch.sizedeg ./ 2;
mR.ppd     = scr.ppdX ./ 2; %because of propixx cutting the screen into 4
mR.degSmoo = 0.2; %in dva
mmR        = JY_VisExptTools('make_smooth_circular_mask', mR);
RIFT.patch.img   = ones( RIFT.patch.sizepix ) .* scr.gray;
RIFT.patch.aper  = 1 - mmR;
RIFT.patch.img0  = RIFT.patch.img .* (1 - mmR) + scr.gray; %aperture = 1-mmR;
RIFT.patch.outer = mmR;


% === After dividing the 1920x1080 frame into four 960x540 frames ===
% -------- Fixation (as a dot) ---------
stim.Fix.size    = round( stim.Fix.size ./ 2 );
stim.Fix.rectFix = cell(4, 1);
for k = 1:4
    stim.Fix.rectFix{k} = CenterRectOnPoint([0 0 stim.Fix.size stim.Fix.size], RIFT.xCtrQuad(k), RIFT.yCtrQuad(k));            
end

% ------- Stimuli of interest ------
% The number of pixels representing the stimulus is halfed
stim.patch.sizepix = ceil( stim.patch.sizepix ./ 2 );
stim.noise         = stim.patch;
% Same logic for the aperture, hence the same operation
m.maskSiz = stim.patch.sizepix;
m.outerR  = stim.patch.sizedeg./2;
m.ppd     = scr.ppdX ./ 2; %because of propixx cutting the screen into 4
m.degSmoo = 0.8; %in dva
mm        = JY_VisExptTools('make_smooth_circular_mask', m);
stim.outerpart = mm; %will only be used when RIFTing luminance
stim.aperture  = 1-mm;

% ------- Dimension Cue ---------
stim.dimcue.ovalrectpix    = stim.dimcue.ovalrectpix ./ 2; %half
stim.dimcue.ovalrects      = cell( 4, 1 );
for k = 1:4 %JY: hard-coded 4
stim.dimcue.ovalrects{k}(:,1) = CenterRectOnPoint(stim.dimcue.ovalrectpix, RIFT.xCtrQuad(k) - 1.75*scr.ppdX/2, RIFT.yCtrQuad(k)); %JY: hard-coded "/2"
stim.dimcue.ovalrects{k}(:,2) = CenterRectOnPoint(stim.dimcue.ovalrectpix, RIFT.xCtrQuad(k) + 1.75*scr.ppdX/2, RIFT.yCtrQuad(k));
end

% --------- Spatial Cue ------------
% % tmpCue4Valid = {'R','D'}; %rectangle vs. diamond
% % stim.spacue.cue4valid = tmpCue4Valid{ mod( str2double(SubjectInfo.SubjectID),2) +1 };
stim.spacue.rtg_rectR= cell(4, 1);
stim.spacue.rtg_rectL= cell(4, 1);
stim.spacue.dm_rectR = cell(4, 1);
stim.spacue.dm_rectL = cell(4, 1);
for k = 1:4 %JY: hard-coded 4
    
    stim.spacue.rtg_rectL{k} = stim.dimcue.ovalrects{k}(:,1);
    stim.spacue.rtg_rectR{k} = stim.dimcue.ovalrects{k}(:,2);
    [xL,~]  = RectCenter(stim.spacue.rtg_rectL{k});
    [xR,~]  = RectCenter(stim.spacue.rtg_rectR{k});
    radius  = floor( stim.dimcue.ovalrectpix(3) ./ 2 .* sqrt(2) );
    
    stim.spacue.dm_rectR{k}(:,1) = [xR, xR-radius, xR, xR+radius];%x-coordinates
    stim.spacue.dm_rectR{k}(:,2) = [RIFT.yCtrQuad(k)+radius, RIFT.yCtrQuad(k), RIFT.yCtrQuad(k)-radius, RIFT.yCtrQuad(k)]; %y-coordinates
    stim.spacue.dm_rectL{k}(:,1) = [xL, xL-radius, xL, xL+radius];%x-coordinates
    stim.spacue.dm_rectL{k}(:,2) = [RIFT.yCtrQuad(k)+radius, RIFT.yCtrQuad(k), RIFT.yCtrQuad(k)-radius, RIFT.yCtrQuad(k)]; %y-coordinates
end


% === Define the triggers ===
global triggers
triggers          = [];
triggers.trlstart = 100;
triggers.dimcue   = 2;
triggers.ICI      = 4;
triggers.spacue   = 8;
triggers.ISI1     = 16;
triggers.S1       = 32;
triggers.ISI2     = 64;
triggers.S2       = 128;


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
texDio   = Screen('MakeTexture', wPtr, RIFT.patch.img0);
reminder = ['Blue for SAME. Red for DIFFERENT.\n\n',stim.spacue.cue4valid ];
for k = 1:4
    drawTextRIFT(wPtr, reminder, k, [], [], [] );
    Screen('DrawTextures', wPtr, [texDio,texDio], [], [RIFT.patch.locRectsR{k};RIFT.patch.locRectsL{k}]');
end
Screen('Flip', wPtr); pause(2); KbWait(ctl.DeviceNum);



%% =========================================================== 
%% ===            calibrate the eye tracker                ===
%% ===========================================================
if ET.use
    
    % Change the sequencer back to the normal non-RIFT mode
    Datapixx('SetPropixxDlpSequenceProgram', 0); %Revert to standard 120Hz refresh rate
    Datapixx('RegWrRd');
    
    % Initialize eye tracker
    [ET.el exitFlag] = rd_eyeLink('eyestart_customize', wPtr, ET.file);
    if exitFlag, return, end
    
    % Calibrate the eye tracker
    [cal exitFlag] = rd_eyeLink('calibrate_customize',wPtr, ET.el);
    if exitFlag, return, end
    
    if EnvInfo.propixx
        switch freqRIFT
            case 1440
                Datapixx('SetPropixxDlpSequenceProgram', 5);
                disp('RIFTing at 1440 Hz!');
            case 480
                Datapixx('SetPropixxDlpSequenceProgram', 2);
                disp('RIFTing at 480 Hz!');
        end
        % Synchronize DATAPixx registers to local register cache
        Datapixx('RegWrRd');
    end
end



%% =========================================================== 
%% ===                 Run the experiment                  ===
%% =========================================================== 
for b = 1:proc.nBlocks
    
    results = cell( proc.nTrialsPerBlock, 1 );
    
    
    % === show block instruction ===
    Screen('FillRect', wPtr, [repmat(scr.gray,[1,3]),1], scr.rect );
    bStartText = ['Waiting for the experimenter to start Block ', num2str(b), '.'];
    for k=1:4, drawTextRIFT(wPtr, bStartText, k, [], [], [] ); end
    Screen('Flip', wPtr); pause(2); KbWait(ctl.DeviceNum);
    
    
    % === start eye-tracking ===
    if ET.use, Eyelink('StartRecording'); end
    
    
    % === buffer ===
    for k = 1:4
        Screen('FrameOval', wPtr, scr.black, stim.Fix.rectFix{k}, 2, 2); %fixation dot
    end
    bStart = Screen('Flip', wPtr);
    
    
    % === loop through the trials ===
    for ii = 1:proc.nTrialsPerBlock
        
        % Define the trial
        trialcfg = trialdef_PredOL_MEG( wPtr, b, ii ); 
        if ii==1
            trialcfg.tStart = bStart + 1 + trialcfg.ITI;
        else
            trialcfg.tStart = results{ii-1}.resp.tEnd + trialcfg.ITI;
        end
        
        % Save the trial configuration parameters
        results{ii}.cfg  = trialcfg;
        results{ii}.resp = [];
        
        % Show trial
        results{ii}.resp = trialfun_PredOL_MEG(wPtr, trialcfg);
        
        % End trial by removing everything but fixation
        for k = 1:4
            Screen('FrameOval', wPtr, scr.black, stim.Fix.rectFix{k}, 2, 2); %fixation dot
        end
        Screen('Flip', wPtr);
        
        % Clear the textures
        Screen('Close',[ results{ii}.cfg.txDioL{:}, results{ii}.cfg.txDioR{:} ]);
        Screen('Close',[ results{ii}.cfg.texL{:}, results{ii}.cfg.texR{:} ]);
    end
    
    % === stop eye-tracking ===
    if ET.use, Eyelink('StopRecording'); end
    
    
    % === end of the block instruction ===
    vec_corr = cellfun(@(x) x.resp.correct, results, 'uniformoutput',true);
    bAcc     = sum(vec_corr) / proc.nTrialsPerBlock;
    bEndText = ['End of Block ', num2str(b), '.\n\n',...
                'You responded correctly on ', num2str(bAcc*100), '% of the trials!\n\n',...
                'Take a break while we save your data...'];
    Screen('FillRect', wPtr, [repmat(scr.gray,[1,3]),1], scr.rect );  
    for k=1:4, drawTextRIFT(wPtr, bEndText, k, [], [], [] ); end
    Screen('Flip', wPtr);
    
    
    % === save data ===
    bFilename = sprintf('b%s_%s', num2str(b), SubjectInfo.datafile);
    save( bFilename, 'results' );
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