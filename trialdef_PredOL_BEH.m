function trialcfg = trialdef_PredOL_BEH( wPtr, b, ii )
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

global stim proc

nCount          = (b-1)* proc.nTrialsPerBlock + ii;
trialcfg        = [];
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

trialcfg.S1 = genBandpassGrating( stim.patch, trialcfg.S1Ori ); 
trialcfg.S2 = genBandpassGrating( stim.patch, trialcfg.S2Ori ); 

trialcfg.M1 = genBandpassPlaid( stim.patch ); %genBandpassPlaidRaw( stim.patch ); 
trialcfg.M2 = genBandpassPlaid( stim.patch ); %genBandpassPlaidRaw( stim.patch ); 

% ========== Prepare textures for drawing ========
% pre-allocate memory
texL = cell( 2, 1 );
texR = cell( 2, 1 );

% input stream on both sides
switch trialcfg.S1Loc
    case 'L', inputL = 'S1'; inputR = 'M1';
    case 'R', inputR = 'S1'; inputL = 'M1';
end
switch trialcfg.S2Loc
    case 'L', inputL = {inputL, 'S2'}; inputR = {inputR, 'M2'};
    case 'R', inputR = {inputR, 'S2'}; inputL = {inputL, 'M2'};
end

% prepare textures for S1 and S2
tmp = stim.patch.patchcon.*stim.aperture;
for ii = 1:2
    inLeft   = eval(['trialcfg.',inputL{ii},'.img']);
    inRight  = eval(['trialcfg.',inputR{ii},'.img']);
    texL{ii} = Screen('MakeTexture', wPtr, inLeft.*tmp+stim.patch.patchlum);
    texR{ii} = Screen('MakeTexture', wPtr, inRight.*tmp+stim.patch.patchlum);
end

% attach the pre-defined textures to the output structure
trialcfg.texL = texL;
trialcfg.texR = texR;


end