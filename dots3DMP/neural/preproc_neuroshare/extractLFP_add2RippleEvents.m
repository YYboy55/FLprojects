function extractLFP_add2RippleEvents(dateRange, varargin)
% extractLFPfromNS2 - Extracts LFP data from .ns2 files for sessions in dateRange, using getData_NS().
%                       Inputs to getData() are supplied by loaded nsEvents file
%                       nsEvents are created by processTrellisData (should be run at the end of every recording experiment)
%                       Saves each resulting extraction alongside the corresponding RippleEvents file.
%
% Required Input:
%   dateRange         - [startDate; endDate] in yyyymmdd format (string or numeric)
%
% Optional Name-Value Inputs:
%   'overwriteFlag'     - (default: false) Overwrite existing RippleEvents file
%   'remoteSessPath'    - (default: hardcoded example) Path to ns2 session directories
%   'localSessInfoPath' - (default: same as remote) Path to local RippleEvents .mat files
%   'useLocal'          - (default: true) Use local RippleEvents if available
%
% current assumptions are:
%     - local nsEvents files for all sessions are in a shared single folder, simply check directly for targetDate file
%     - .ns2 and nsEvents files are in a session subfolder together, named after the session date, i.e. enter to targetDate subfolder then grab files
%     - nsEvents file names contain the targetDate and end in 'RippleEvents.mat'
%     - multiple ns2 and nsEvents files can exist for the same targetDate, these are from same day but different recording sets and are thus treated separately (read processTrellisData for more details)



addpath(genpath('C:\Users\yhaile2\Documents\AcademicRelated\CODE_Projects\GitHubCodes\Fetschlab\FLprojects\dots3DMP\neural')) % make sure getData on path
addpath('C:\Users\yhaile2\Documents\AcademicRelated\CODE_Projects\MATLAB\NeuroShareFunctions'); % neuroShare functions

% ---> ACTIVATE LIBRARY <---    once per matlab instance to ensure neuroshare applications are functional! 
% status = ns_SetLibrary('C:\Users\yhaile2\Documents\AcademicRelated\CODE_Projects\MATLAB\NeuroShareFunctions\nsNEVLibrary.dll') % 0 means library succesfully loaded/initialized

% Input parsing
p = inputParser;
p.addRequired('dateRange');
p.addParameter('overwriteFlag', false, @islogical);
p.addParameter('remoteSessPath', 'Z:\fetschlab\data\lucio\lucio_neuro', @ischar);
p.addParameter('localSessInfoPath', 'C:\Users\yhaile2\Documents\AcademicRelated\CODE_Projects\Data\Lucio\Neural\SessionInfoFiles_from_GetData', @ischar);
p.addParameter('useLocal', false, @islogical);
p.parse(dateRange, varargin{:});

overwriteFlag = p.Results.overwriteFlag;
remoteSessPath = p.Results.remoteSessPath;
localSessInfoPath = p.Results.localSessInfoPath;
useLocal = p.Results.useLocal;

% Convert dateRange to string
if isnumeric(dateRange), dateRange = string(dateRange); end

try
    startDate = min(datetime(dateRange, 'InputFormat', 'yyyyMMdd'));
    endDate   = max(datetime(dateRange, 'InputFormat', 'yyyyMMdd'));
catch

    error('Input date format not supported, must be yyyyMMdd')
end

% select valid session folders within specified date range
sessFolders = dir(remoteSessPath);
sessFolders = sessFolders([sessFolders.isdir]);
sessFolders = sessFolders(~ismember({sessFolders.name}, {'.','..'}));
% Parse folder names as dates
allSessDates = datetime({sessFolders.name}, 'InputFormat', 'yyyyMMdd', 'Format', 'yyyyMMdd', 'Locale', 'en_US');
% Filter dates by set dateRange
validIdx = allSessDates >= startDate & allSessDates <= endDate;
filteredDates = allSessDates(validIdx);

for sess = 1:length(filteredDates)
    targetDate = char(filteredDates(sess), 'yyyyMMdd');
    currSessPath = fullfile(remoteSessPath,targetDate); % ns2 always from server (i.e. use remote path)

    % Locate RippleEvents files
    if useLocal
        nsEvFiles_local = dir(fullfile(localSessInfoPath, sprintf('*%s*_RippleEvents.mat', targetDate)));
        if ~isempty(nsEvFiles_local)
            nsEvFiles = {nsEvFiles_local.name};
            nsEvPath = localSessInfoPath;
        else
            nsEvFiles = {dir(fullfile(currSessPath, sprintf('*%s*_RippleEvents.mat', targetDate))).name}; % this replaces below, jsut make sure it works as expected
%             nsEvFiles = dir(fullfile(currSessPath,sprintf('*%s*_RippleEvents.mat',targetDate)));
%             nsEvFiles = {nsEvFiles.name};
            nsEvPath = currSessPath;
        end
    else
        nsEvFiles = {dir(fullfile(currSessPath, sprintf('*%s*_RippleEvents.mat', targetDate))).name};
        nsEvPath = currSessPath;
    end

    ns2Files = {dir(fullfile(currSessPath, '*.ns2')).name}; % same deal here
%     ns2Files = dir(fullfile(currSessPath,'*.ns2'));
%     ns2Files =  {ns2Files.name};

% load files for each rec file in curr sess i.e. each 'set'
    for set = 1:length(ns2Files)
        S = load( fullfile(nsEvPath,nsEvFiles{set}) );
        if ~isfield(S, 'nsEvents')
            warning('Missing nsEvents for %s', [targetDate '_00' num2str(set)]);
            continue;
        end
        nsEvents = S.nsEvents;
        % getData inputs
        channels = str2double(nsEvents.analogInfo.label);
        ns2FilePath = fullfile(currSessPath, ns2Files{set});
        % init output storage
        nCh = length(channels);
        lfpDataCell = cell(nCh, 1);
        chanLabels  = cell(nCh, 1);
        % temporary cell to store fs and time
        tempData = cell(nCh, 1);

        % to avoid loading same file for each ch load once here, and pass hFile to getData
%         [ns_status, hFile] = ns_OpenFileSJ(ns2FilePath, 'single');
%         if ns_status ~= 0
%             error('Failed to open ns2 file.');
%         end
%         % in loop over ch
%         nsData = getdata_NS(hFile, 'LFP', ch);

        for c = 1:length(channels)
            ch = channels(c);

            nsData = getdata_NS(ns2FilePath, 'LFP', ch);

            lfpDataCell{c} = nsData.analogData;
            chanLabels{c}  = nsData.analogInfo.EntityLabel; % same as input nsEvents chLabel, no?

            tempData{c} = struct('fs', nsData.analogInfo.SampleRate, ...
                'time', nsData.analogTime);
        end

        % same for all ch, just need one
        fs = tempData{1}.fs;
        lfpTime = tempData{1}.time;

        % Convert cell array to matrix
        nSamples = length(lfpDataCell{1});
        lfpData = zeros(nCh, nSamples);
        for c = 1:nCh
            lfpData(c, :) = lfpDataCell{c};
        end

        % Store in struct
        lfpStruct = struct('date', targetDate, 'set', set, 'data', lfpData, ...
            'time', lfpTime, 'fs', fs, 'chanLabels', {chanLabels}, 'chanIDs', channels);

        % for now seems fine to store LFP as stand alone variable in same file. Make sure downstream scripts load all vars from RippleEvents
% %         nsEvents.LFP = lfpStruct; % do we want as a new field or alternate way of storing?

        currFile = fullfile(nsEvPath, nsEvFiles{set});
        [~, baseName, ~] = fileparts(currFile);

        if overwriteFlag
            saveFile = fullfile(currSessPath, baseName);  %#ok<UNRCH>   % Overwrite the original RippleEvents file
        else
            saveFile = fullfile(currSessPath, [baseName '_wLFP.mat']);  % Save new file, adding '_wLFP' suffix
        end

        fprintf('Saving LFP data to: %s\n', saveFile);
        save(saveFile, 'nsEvents', 'lfpStruct', '-v7.3');

        clear lfpStruct
    end
end

