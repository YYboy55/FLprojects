% clean monkey data structure
% now incorporated components missing from predecessor 'dots3DMP_cleanMonkeyData' % YYY, 2025-05-29 %

% struct data has fields:
% filename
% subj: subject code
% choice: 1=left, 2=right, nan = fixation break or otherwise invalid
% heading: angle in deg relative to straight ahead, positive = rightward
% coherence: dot motion coherence aka visual reliability
% modality: stimulus modality: 1=ves, 2=vis, 3=comb
% delta: conflict angle in deg, positive means visual right, vestib left
% correct: was the choice correct, 1 or 0
% conf: confidence rating via saccadic end-point to color-bar target
%       (or in other datasets, PDW)

clear all; close all
conftask = 2; % 2 == binary PDW
RTtask = 1; % 1==true
today    = str2double(datestr(now,'yyyymmdd')); % datetime("today","InputFormat","uuuu-MM-dd")


paradigm = 'dots3DMP';
% dateRange = 20210315:20210805; % RT
% dateRange = 20211101:20220809; % RT
% dateRange = 20230501:20250616; %20220512:20230605; %Start and end date of sessions in file of interest
dateRange = 20220512:20220514; % ADC eyePos data
%% Load data file
subject = 'lucio';
folder = 'C:\Users\yhaile2\Documents\AcademicRelated\CODE_Projects\Data\Lucio\Behavior\ConcatBehavSessFiles';
store_folder = 'C:\Users\yhaile2\Documents\AcademicRelated\CODE_Projects\Data\Lucio\Behavior\ConcatBehavSessFiles\CleanedDataFiles';

% subject = 'zarya';
% folder = 'C:\Users\yhaile2\Documents\AcademicRelated\CODE_Projects\Data\Zarya\Behavior\ConcatBehavSessFiles';
% store_folder = 'C:\Users\yhaile2\Documents\AcademicRelated\CODE_Projects\Data\Zarya\Behavior\ConcatBehavSessFiles\CleanedDataFiles';

% file = [subject '_' num2str(dateRange(1)) '-' num2str(dateRange(end)) '.mat'];
file = [subject '_' num2str(dateRange(1)) '-' num2str(dateRange(end))];
loadFile = fullfile(folder,file);
load(loadFile,'data');

try data = rmfield(data,'amountRewardLowConf'); catch, end
try data = rmfield(data,'amountRewardHighConf'); catch, end

% new 04/2020: selecting 'good' data (esp RT)

%% some new useful vars
for k = 1:length(data.filename)
    data.date(k,1) = str2double(data.filename{k}(6:13)); % 6:13 == selecting only the date of the file name
    data.subjDate{k,:} = [data.subj{k} data.filename{k}(6:13)]; % Variable now lists filenames as zaryaYYYYMMdd
end


%% Some manual excludes e.g. of bad sessions, non-RT/PDW

excludes_filename = {};
excludes_date = [];

% remove fixation breaks (indicated by nans), or manually excluded filenames
removethese = isnan(data.choice) | isnan(data.PDW) ; % | isnan(data.RT) | isinf(data.RT) ;
removethese = removethese | ismember(data.filename,excludes_filename) | ismember(data.date,excludes_date);

fnames = fieldnames(data);
for F = 1:length(fnames)
    if strcmp(fnames{F},'velProfile'), continue; end % since not saved every trial indexing is off
    if strcmp(fnames{F},'posTraj'), continue; end % since not saved every trial indexing is off
    eval(['data.' fnames{F} '(removethese) = [];']);
end

% should do the trial number based exclusion here, once brfixes are removed

% quick look at blocks, for when some need to be excluded
blocks = unique(data.filename); % Trials in same block have same file name
nTrialsByBlock = nan(length(blocks),1);
for u = 1:length(blocks)
    nTrialsByBlock(u) = sum(ismember(data.filename,blocks(u)));
end

% discard blocks with <N good trials
N = 50;
removethese = ismember(data.filename,blocks(nTrialsByBlock<N));

for F = 1:length(fnames)
    fieldData = data.(fnames{F});
    if ~isempty(fieldData) && numel(removethese) == numel(fieldData) %  skip empty fields and fields with ntrials < all trials (aka velProfile currently)
        data.(fnames{F})(removethese) = [];
    end
end

% quick look at blocks again
blocks = unique(data.filename);
nTrialsByBlock = nan(length(blocks),1);
for u = 1:length(blocks)
    nTrialsByBlock(u) = sum(ismember(data.filename,blocks(u)));
end


%% cull data

mods = unique(data.modality);

%  consolidate deltas
data.delta(data.delta<0) = -3;
data.delta(data.delta>0) = 3;
deltas = unique(data.delta); % aka conflict angle

% simplify cohs (collapse similar ones)
data.coherence(data.coherence<0.5) = 0.3;
data.coherence(data.coherence>0.5) = 0.7;
cohs = unique(data.coherence);

% the coh assigned to vestib trials (as a placeholder) depends on which
% cohs were shown in a particular block, so we need to standardize it:
data.coherence(data.modality==1) = cohs(1);

data.heading(abs(data.heading)<0.01) = 0;
hdgs = unique(data.heading);
% consolidate heading values
hdg_ranges = [0 1 2 3 4 5 8 10 12];
newhdgvals = [0 1.5 1.5 3 3 6 6 12 12]; % match length of hdg_ranges
hdg_index = interp1(hdg_ranges, 1:numel(hdg_ranges), abs(data.heading), 'nearest');
data.heading_recoded = newhdgvals(hdg_index)' .* sign(data.heading);

% remove the rest
removethese = ~ismember(data.heading,hdgs) | ~ismember(data.coherence,cohs);
for F = 1:length(fnames)
    eval(['data.' fnames{F} '(removethese) = [];']);
end

% Shouldn't be an issue for current datasets
% fix data.correct (see function description for issue!)
% data.correct = dots3DMPCorrectTrials(data.choice,data.heading,data.delta);

% remove one target trials
% removethese = data.oneTargChoice | data.oneTargConf;
% for F = 1:length(fnames)
%     eval(['data.' fnames{F} '(removethese) = [];']);
% end

% try data = rmfield(data,'reward'); catch, end
% try data = rmfield(data,'subj'); catch, end
% try data = rmfield(data,'oneTargChoice'); catch, end
% try data = rmfield(data,'oneTargConf'); catch, end
try data = rmfield(data,'TargMissed'); catch, end
try data = rmfield(data,'subjDate'); catch, end
try data = rmfield(data,'insertTrial'); catch, end
try data = rmfield(data,'confRT'); catch, end

%% check block and day trial counts
blocks = unique(data.filename);
nTrialsByBlock = nan(length(blocks),1);
for u = 1:length(blocks)
    nTrialsByBlock(u) = sum(ismember(data.filename,blocks(u)));
end

dates = unique(data.date);
nTrialsByDate = nan(length(dates),1);
for u = 1:length(dates)
    nTrialsByDate(u) = sum(ismember(data.date,dates(u)));
end

%% save it
save(fullfile(store_folder,[file '_clean.mat']),'data','nTrialsByDate','nTrialsByBlock');
% save(fullfile(store_folder,[file '_clean.mat']),'data'); % store just 'data'
% save(fullfile(store_folder,[file(1:end-4) '_clean.mat']),'data','nTrialsByDate','nTrialsByBlock'); % remove '.mat' if included as part of file name set at top
