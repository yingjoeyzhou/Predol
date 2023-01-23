function trialcfg = trialdef_PredOL_MEG( wPtr, b, ii )
% Define (and prepare) for the trial.
%
% INPUT:
%   wPtr   : window pointer.
%   b      : block index.
%   ii     : trial (within-the-block) index.
%
% OUTPUT:
%   trialcfg: structure, fields contain trial-specific info.
%
% Joey Zhou (Nov 2022)

global stim proc RIFT scr

nCount          = (b-1)* proc.nTrialsPerBlock + ii;
trialcfg        = [];
% trialcfg.ISI2   = stim.dur.ISI2( proc.idxISI(nCount) );%randi( stim.dur.ISI2.*1000 ) / 1000;
trialcfg.ITI    = randi( stim.dur.ITI.*1000 ) / 1000; %for programing convenience, ITI is the preiod *before* trialstart

trialcfg.task   = proc.Task(nCount);
trialcfg.cue    = proc.CueVal(nCount);
trialcfg.S1Loc  = proc.S1Loc(nCount);
trialcfg.S2Loc  = proc.S2Loc(nCount);
if strcmpi( trialcfg.cue, 'V' )
    trialcfg.cue = trialcfg.S1Loc;
end

trialcfg.S1Ori = proc.S1Ori(nCount);
trialcfg.S2Ori = proc.S2Ori(nCount);

trialcfg.S1 = genBandpassGratingRaw( stim.patch, trialcfg.S1Ori );
trialcfg.S2 = genBandpassGratingRaw( stim.patch, trialcfg.S2Ori );

trialcfg.M1 = genBandpassPlaidRaw( stim.patch ); %genBandpassFilteredNoiseRaw( stim.noise );
trialcfg.M2 = genBandpassPlaidRaw( stim.patch ); %genBandpassFilteredNoiseRaw( stim.noise );

% ======== OLD ========
% trialcfg.bMask = genBandpassFilteredNoiseRaw( stim.noise );
% =====================
% % % % trialcfg.bMask     = trialcfg.S1; %arbitrary
% % % % trialcfg.bMask.img = trialcfg.bMask.img.*0;
trialcfg.bMask     = trialcfg.S1; %arbitrary
trialcfg.bMask.img = ones( size(trialcfg.S1.img) );

% Define the frequency-tagging signal
switch RIFT.whichsig(nCount)
    case 1 %signals from s1 matrix is applied to the left side
        trialcfg.RIFT_L = RIFT.s1(:, RIFT.indexsig(nCount));
        trialcfg.RIFT_R = RIFT.s2(:, RIFT.indexsig(nCount));
    case 2
        trialcfg.RIFT_R = RIFT.s1(:, RIFT.indexsig(nCount));
        trialcfg.RIFT_L = RIFT.s2(:, RIFT.indexsig(nCount));
end


% ============= debug
%{
dur = 1.6;
fs  = 1440;
t   = 0:(1/fs):(dur-1/fs);
trialcfg.RIFT_L = defineBroadbandLuminance( t, 80 );
trialcfg.RIFT_R = defineBroadbandLuminance( t, 70 );
%}
% ============= debug


% ========== Prepare textures for drawing ========
% reshape the freq-tagged signal (i.e., contrast)
tagL = reshape( trialcfg.RIFT_L, [12, numel(trialcfg.RIFT_L)./12] ) .*0.5 + 0.25;
tagR = reshape( trialcfg.RIFT_R, [12, numel(trialcfg.RIFT_R)./12] ) .*0.5 + 0.25;

disp(range(tagL(:)) );

% adjust the range: tag luminance
trialcfg.S1.img = trialcfg.S1.img .* 0.5;
trialcfg.S2.img = trialcfg.S2.img .* 0.5;
trialcfg.M1.img = trialcfg.M1.img .* 0.5;
trialcfg.M2.img = trialcfg.M2.img .* 0.5;

% pre-allocate memory
texL = cell( RIFT.nFrameTot, 1 );
texR = cell( RIFT.nFrameTot, 1 );
txDioL = cell( RIFT.nFrameTot, 1 );
txDioR = cell( RIFT.nFrameTot, 1 );

% input stream on both sides
switch trialcfg.S1Loc
    case 'L', inputL = 'S1'; inputR = 'M1';
    case 'R', inputR = 'S1'; inputL = 'M1';
end
switch trialcfg.S2Loc
    case 'L', inputL = {inputL, 'S2'}; inputR = {inputR,'M2'};
    case 'R', inputR = {inputR, 'S2'}; inputL = {inputL,'M2'};
end

cnt = 0;

for iFlip = 1:RIFT.nFrameTot
    tmpL = tagL(:,iFlip);
    tmpR = tagR(:,iFlip);
    tmpL = reshape( tmpL, [1, 4, 3] ); %from 1st to 3rd column, Reds -> Greens -> Blues;
    tmpR = reshape( tmpR, [1, 4, 3] ); %from 1st to 3rd column, Reds -> Greens -> Blues;
    
    if iFlip <= RIFT.nFrame0
        inLeft  = trialcfg.bMask.img;
        inRight = trialcfg.bMask.img;
        
        stimOn = 0;
        
    elseif (iFlip>RIFT.nFrame0) & (iFlip<=(RIFT.nFrame0+RIFT.nFrameS1))
        inLeft = eval(['trialcfg.',inputL{1},'.img']);
        inRight= eval(['trialcfg.',inputR{1},'.img']);
        
        stimOn = 1;
        
    elseif (iFlip>(RIFT.nFrame0+RIFT.nFrameS1)) & (iFlip<=(RIFT.nFrame0+RIFT.nFrameS1+RIFT.nFrameISI))
        inLeft  = trialcfg.bMask.img;
        inRight = trialcfg.bMask.img;
        
        stimOn = 0;
        
    elseif (iFlip>(RIFT.nFrame0+RIFT.nFrameS1+RIFT.nFrameISI))
        inLeft = eval(['trialcfg.',inputL{2},'.img']);
        inRight= eval(['trialcfg.',inputR{2},'.img']);
        
        stimOn = 1;
    end
    
    for k = 1:4 %The four quardrants
        
        cnt = cnt + 1;
        
        switch stimOn
            case 0
                imgL = ( inLeft .* tmpL(1,k,:)).* stim.aperture + stim.outerpart.*stim.patch.patchlum;
                imgR = ( inRight .* tmpR(1,k,:)).* stim.aperture + stim.outerpart.*stim.patch.patchlum;
            case 1
                imgL = ( inLeft + tmpL(1,k,:)).* stim.aperture + stim.outerpart.*stim.patch.patchlum;
                imgR = ( inRight + tmpR(1,k,:)).* stim.aperture + stim.outerpart.*stim.patch.patchlum;
        end
        
        dioL = ( RIFT.patch.img + tmpL(1,k,:)-0.25 ).* RIFT.patch.aper + RIFT.patch.outer.*stim.patch.patchlum;
        dioR = ( RIFT.patch.img + tmpR(1,k,:)-0.25 ).* RIFT.patch.aper + RIFT.patch.outer.*stim.patch.patchlum;
        
        % disp( [max(imgL(:)), min(imgL(:))] );
        
        texL{cnt} = Screen('MakeTexture',wPtr,imgL);
        texR{cnt} = Screen('MakeTexture',wPtr,imgR);
        txDioL{cnt} = Screen('MakeTexture',wPtr,dioL);
        txDioR{cnt} = Screen('MakeTexture',wPtr,dioR);
        
    end
end

trialcfg.texL = texL;
trialcfg.texR = texR;
trialcfg.txDioL = txDioL;
trialcfg.txDioR = txDioR;


end