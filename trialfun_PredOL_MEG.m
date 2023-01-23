function resp = trialfun_PredOL_MEG( wPtr, trialcfg )

global devMode scr stim ctl RIFT triggers ET

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
drawDimCue_RIFT( wPtr, 'P' ); %P for placeholder
% flip to screen
tStart = Screen('Flip',wPtr,trialcfg.tStart - scr.ifi/2);
if ET.use, Eyelink('Message', 'EVENT_TSTART'); end
if devMode==false, sendtrig(triggers.trlstart); end


% ------- trial start with dimension cue -------
% draw dimension cue
drawDimCue_RIFT( wPtr, dimcue );
% flip to screen
tDimCue = Screen('Flip', wPtr, tStart + stim.dur.precue - scr.ifi/2);
if ET.use, Eyelink('Message','EVENT_DIMCUE'); end
if devMode==false, sendtrig(triggers.dimcue); end


% --------- cue-to-cue ISI -----------
% draw dimension cue
drawDimCue_RIFT( wPtr, 'P' ); %P for placeholder
% flip to screen
tISI = Screen('Flip', wPtr, tDimCue + stim.dur.dimcue - scr.ifi/2);
if devMode==false, sendtrig(triggers.ICI); end


% --------- spatial cue ----------
% draw spatial cue
drawSpatialCue_RIFT(wPtr, trialcfg.cue);
% flip to screen
tSpaCue = Screen('Flip', wPtr, tISI + stim.dur.ISI - scr.ifi/2);
if ET.use, Eyelink('Message','EVENT_SPACUE'); end
if devMode==false, sendtrig(triggers.spacue); end


% ------------ start RIFTing -----------
% The cue-to-stim ISI, S1, stim-to-stim ISI, and S2 are RIFTing.
frame_cnt       = 0;
event_cnt       = 0;
event_onsets    = nan(1,4); %JY: hard-coded
event_triggers  = [triggers.ISI1, triggers.S1, triggers.ISI2, triggers.S2];
event_names     = {'EVENT_ISI1', 'EVENT_S1', 'EVENT_ISI2', 'EVENT_S2'};
event_prevdur   = [stim.dur.spacue, stim.dur.cue2stim, stim.dur.s1, stim.dur.ISI];

prev_onset = tSpaCue;

fliptime_RIFT = nan(1, RIFT.nFrameTot);

for iF = 1:RIFT.nFrameTot

% % %     % for debugging
% % %     disp(iF);
    
    % draw stimuli
    for k = 1:4 %the four quardrants
        frame_cnt = frame_cnt + 1;
        Screen('DrawTextures', wPtr, [trialcfg.texL{frame_cnt}, trialcfg.texR{frame_cnt}, trialcfg.txDioR{frame_cnt}, trialcfg.txDioL{frame_cnt}],[],...
                [RIFT.rect_L{k}; RIFT.rect_R{k}; RIFT.patch.locRectsR{k}; RIFT.patch.locRectsL{k}]',[],[]);
        Screen('FrameOval', wPtr, scr.black, stim.Fix.rectFix{k}, 2, 2); %fixation dot
        Screen('FrameOval', wPtr, [0,0,0; 0,0,0]', [RIFT.rect_L{k}; RIFT.rect_R{k}]', 2, 2); %contour of the stimuli
    end

    % send trigger if the current frame is the first one of an event
    if any(RIFT.frameOnset==iF) %when this is the first frame of an event
        event_cnt = event_cnt + 1;
        event_vbl = Screen('Flip',wPtr, prev_onset+event_prevdur(event_cnt)-scr.ifi/2);
        
        prev_onset = event_vbl;
        
        event_onsets(event_cnt) = event_vbl;
        
        if devMode==false, sendtrig( event_triggers(event_cnt) ); end
        if ET.use, Eyelink('Message',event_names{event_cnt}); end
        
        fliptime_RIFT(iF) = event_vbl;
        
% % %         % for debugging
% % %         disp( event_names{event_cnt} );
        
    else %when this is NOT the first frame of an event
        vbl = Screen('Flip',wPtr);
        fliptime_RIFT(iF) = vbl;
        
        % pause(0.1)
    end
end


% ------------ response ------------
for k = 1:4 %The four quardrants
    Screen('FrameOval', wPtr, scr.black, stim.Fix.rectFix{k}, 2, 2); %fixation dot
    Screen('FrameOval', wPtr, [0,0,0; 0,0,0]', [RIFT.rect_L{k}; RIFT.rect_R{k}]', 2, 2); %contour
end
% flip to the screen
tStimOff = Screen('Flip',wPtr, prev_onset + stim.dur.s2 - scr.ifi/2); %JY: stim.dur.mask was set to zero though
% check response
if resp.respCheck == 0
    r          = resp;
    r.timestop = tStimOff + stim.dur.MaxResp - scr.ifi/2;
    resp = JY_VisExptTools('get_keyboard_response', r);
end

% --------- output ------------
resp.tEvents  = [tStart, tDimCue, tISI, tSpaCue, event_onsets, tStimOff];
resp.tEnd     = min( [r.timestop, resp.respTime] );
resp.flipRIFT = fliptime_RIFT;


end


