function resp = trialfun_PredOL_BEH( wPtr, trialcfg )

global devMode scr stim ctl ET

% output structure
resp           = [];
resp.respCheck = 0;

% define the dimension cue text
switch trialcfg.task
    case 'S' 
        dimcue = 'Location'; 
        if strcmpi( trialcfg.S1Loc, trialcfg.S2Loc )
            resp.answer = ctl.keySame;
        else
            resp.answer = ctl.keyDiff;
        end
    case 'O' 
        dimcue = 'Orientation';
        if strcmpi( trialcfg.S1Ori, trialcfg.S2Ori )
            resp.answer = ctl.keySame;
        else
            resp.answer = ctl.keyDiff;
        end
end

% ---------- trial start with cue placeholders -----------
% draw dimension cue
drawDimCue( wPtr, 'P' ); %P for placeholder
% draw stimuli contours
Screen('FrameOval', wPtr, [0,0,0; 0,0,0]', [stim.loc.rect_R; stim.loc.rect_L]');
% flip to screen
tStart = Screen('Flip',wPtr,trialcfg.tStart - scr.ifi/2);
if ET.use, Eyelink('Message', 'EVENT_TSTART'); end


% ------- trial start with dimension cue -------
% draw dimension cue
drawDimCue( wPtr, dimcue );
% draw stimuli contours
Screen('FrameOval', wPtr, [0,0,0; 0,0,0]', [stim.loc.rect_R; stim.loc.rect_L]');
% flip to screen
tDimCue = Screen('Flip', wPtr, tStart + stim.dur.precue - scr.ifi/2);
if ET.use, Eyelink('Message','EVENT_DIMCUE'); end


% --------- cue-to-cue ISI -----------
% draw dimension cue
drawDimCue( wPtr, 'P' ); %P for placeholder
% draw stimuli contours
Screen('FrameOval', wPtr, [0,0,0; 0,0,0]', [stim.loc.rect_R; stim.loc.rect_L]');
% flip to screen
tISI = Screen('Flip', wPtr, tDimCue + stim.dur.dimcue - scr.ifi/2);
if devMode==false, sendtrig(triggers.ICI); end


% --------- spatial cue ----------
% draw spatial cue
drawSpatialCue(wPtr, trialcfg.cue);
% draw stimuli contours
Screen('FrameOval', wPtr, [0,0,0; 0,0,0]', [stim.loc.rect_R; stim.loc.rect_L]');
% flip to screen
tSpaCue = Screen('Flip', wPtr, tISI + stim.dur.ISI - scr.ifi/2);
if ET.use, Eyelink('Message','EVENT_SPACUE'); end


% ------------- ISI -----------------
%fixation dot
Screen('FrameOval', wPtr, scr.black, stim.Fix.rectFix, 4, 4); 
% draw stimuli contours
Screen('FrameOval', wPtr, [0,0,0; 0,0,0]', [stim.loc.rect_R; stim.loc.rect_L]');
tISI1 = Screen('Flip', wPtr, tSpaCue + stim.dur.spacue - scr.ifi/2);


% ------------- S1 ---------------
%fixation dot
Screen('FrameOval', wPtr, scr.black, stim.Fix.rectFix, 4, 4); 
% draw stimuli
Screen('DrawTextures', wPtr, [trialcfg.texR{1}, trialcfg.texL{1}], [],...
                            [stim.loc.rect_R; stim.loc.rect_L]',[],[]);
% draw stimuli contours
Screen('FrameOval', wPtr, [0,0,0; 0,0,0]', [stim.loc.rect_R; stim.loc.rect_L]');
tS1 = Screen('Flip', wPtr, tISI1 + stim.dur.cue2stim - scr.ifi/2);


% ------------ ISI ---------------
%fixation dot
Screen('FrameOval', wPtr, scr.black, stim.Fix.rectFix, 4, 4); 
% draw stimuli contours
Screen('FrameOval', wPtr, [0,0,0; 0,0,0]', [stim.loc.rect_R; stim.loc.rect_L]');
tISI2 = Screen('Flip', wPtr, tS1 + stim.dur.s1 - scr.ifi/2);


% ------------ S2 ----------------
%fixation dot
Screen('FrameOval', wPtr, scr.black, stim.Fix.rectFix, 4, 4); 
% draw stimuli
Screen('DrawTextures', wPtr, [trialcfg.texR{2}, trialcfg.texL{2}], [],...
                            [stim.loc.rect_R; stim.loc.rect_L]',[],[]);
% draw stimuli contours
Screen('FrameOval', wPtr, [0,0,0; 0,0,0]', [stim.loc.rect_R; stim.loc.rect_L]');
tS2 = Screen('Flip', wPtr, tISI2 + stim.dur.ISI - scr.ifi/2);


% ------------ response ------------
%fixation dot
Screen('FrameOval', wPtr, scr.black, stim.Fix.rectFix, 4, 4);
% draw stimuli contours
Screen('FrameOval', wPtr, [0,0,0; 0,0,0]', [stim.loc.rect_R; stim.loc.rect_L]');
tStimOff = Screen('Flip',wPtr, tS2 + stim.dur.s2 - scr.ifi/2); 
% check response
if resp.respCheck == 0
    r          = resp;
    r.timestop = tStimOff + stim.dur.MaxResp - scr.ifi/2;
    resp = JY_VisExptTools('get_keyboard_response', r);
end

% --------- output ------------
resp.tEvents  = [tStart, tDimCue, tISI, tSpaCue, tISI1, tS1, tISI2, tS2, tStimOff];
resp.tEnd     = min( [r.timestop, resp.respTime] );


end