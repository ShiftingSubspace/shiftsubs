function [patterns,markerDict] = patternGen(protocol,op_mode)
% Generate correct patterns (and MarkerDict) for each obj based on (RGM's UTCS/TEMPO) protocol and operation_mode. 
%   Protocol: "UTCS" or "TEMPO"
%   op_mode: "Ex/ObsDelay" or "Ex/Obsnodelay"
if ~exist('op_mode', 'var')
    op_mode = 'ExDelay';
end
if isequal(op_mode, "ExDelay")
    op_code = 34;
elseif isequal(op_mode, "ObsDelay")
    op_code = 35;
elseif isequal(op_mode, "Ex")
    op_code = 32;
end

% todo1: decide which protocol to use based on Monkey & Date. Generally,
% before 2021-04-12 TEMPO, after that UTCS. Monkey I&T have special cases
% (Email from Marc on Nov.29). T_20220526 & I_20210129. (Copy to QuickNotes
% OneNote, ur.sam@outlook.com)
%--> solution: OUTSIDE this function, a new function that specify protocol
%as output. INSIDE this function, ADD more entry options for diff protocol,
%i.e., two special cases. 

numCond = 4;
obj = [8 1 128 64];
patterns = cell(1,numCond);
if isequal(protocol, "UTCS")
    markerDict = struct(...
    'init', Marker(251,'.', '.', '.', '.', '.', '.', '.', '.',31),...
    'inston', 210, ...
    'instoff', 211,...
    'gocue', 245,...
    'mvon', 59 ,...
    'obcontact', 246,...
    'holds', 248,...
    'holde',119 ); % ref: RGM_UTCS_TaskEpochs.pdf (withDelay)

    
    for i=1:numCond   % execution patterns
        head_code = markerDict.init.replace_at(markerDict.init.elements_count-1, obj(i));    
    %     patterns{i} = Pattern(head_code, op_code, markerDict.inston,...
    %         markerDict.mvon, 60, markerDict.obcontact, markerDict.holds, markerDict.holde, 120, 121, 253);
        patterns{i} = Pattern(head_code, op_code, markerDict.inston, markerDict.instoff, markerDict.gocue,...
            markerDict.mvon, 60, markerDict.obcontact, markerDict.holds, markerDict.holde, 120, 121, 253);
    end
elseif isequal(protocol, "TEMPO")
    markerDict = struct(...
    'init', Marker(251, '.', '.', '.', '.', '.', '.', '.', '.', '.','.', 31),...
    'inston', 245, ...
    'instoff', 50,...
    'gocue', 54,...
    'mvon', 64 ,...
    'obcontact', 246,...
    'holds', 248,...
    'holde',70 ); % ref: TempoCodeSequenceRGM.xlsx (withDelay)
    for i=1:numCond
        patterns{i} = Pattern(markerDict.init, op_code, 40, markerDict.inston, obj(i), 128,...
            markerDict.instoff, markerDict.gocue, markerDict.mvon, markerDict.obcontact,...
            markerDict.holds, markerDict.holde, 119, 253);
    end

elseif isequal(protocol, "T_20220526") % prtc for T's special case on 20220526, UTCS-alike, Q: dots diff in init?
    markerDict = struct(...
    'init', Marker(85,'.', '.', '.', '.', '.', '.', '.', '.', '.',31),...
    'inston', 210, ...
    'instoff', 211,...
    'gocue', 245,...
    'mvon', 59 ,...
    'obcontact', 246,...
    'holds', 248,...
    'holde',119 ); % ref: UTCS+email(Nov-29-2022 from Marc) - 85(251), 250(253)

    
    for i=1:numCond   
        head_code = markerDict.init.replace_at(markerDict.init.elements_count-1, obj(i));    
    %     patterns{i} = Pattern(head_code, op_code, markerDict.inston,...
    %         markerDict.mvon, 60, markerDict.obcontact, markerDict.holds, markerDict.holde, 120, 121, 253);
        patterns{i} = Pattern(head_code, op_code, markerDict.inston, markerDict.instoff, markerDict.gocue,...
            markerDict.mvon, 60, markerDict.obcontact, markerDict.holds, markerDict.holde, 120, 121, 250);
    end

elseif isequal(protocol, "I_20210129") 
    % mvon is 60(special, convention is  59)
    % no 246 
%     startTrialDirectionMarkers{p} = Marker(251,'.', '.', '.', '.', '.', '.', '.', targets(p),'.');   
%         %                                                1         2        3   4   5  6  7   8   9  10  11  12  13                                                    
%          patterns{p} = Pattern(startTrialDirectionMarkers{p},ExeObsContext,210,211,245,60,60,248,119,120,121,253);


    markerDict = struct(...
    'init', Marker(251,'.', '.', '.', '.', '.', '.', '.', '.',31),...
    'inston', 210, ...
    'instoff', 211,...
    'gocue', 245,...
    'mvon', 60 ,...
    'obcontact', 246,...
    'holds', 248,...
    'holde',119 ); % ref: UTCS+ email mvon change, remove obcontact

    
    for i=1:numCond   % execution patterns
        head_code = markerDict.init.replace_at(markerDict.init.elements_count-1, obj(i));    
    %     patterns{i} = Pattern(head_code, op_code, markerDict.inston,...
    %         markerDict.mvon, 60, markerDict.obcontact, markerDict.holds, markerDict.holde, 120, 121, 253);
        patterns{i} = Pattern(head_code, op_code, markerDict.inston, markerDict.instoff, markerDict.gocue,...
            markerDict.mvon, 60, markerDict.holds, markerDict.holde, 120, 121, 253);
    end

end

