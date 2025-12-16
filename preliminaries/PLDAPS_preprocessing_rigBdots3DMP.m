% grabs PLDAPS data from NAS as specified, then loads desired variables
% from .PDS files into a simpler data structure

% CF, adapted from offline analysis wrapper born 10-3-18

% VERSION 2.0: 08-30-19
% New toggles: 07-25-2025
%  - include eyePos data (which also automatically includes event time vectors)
%  - include idealized motion trajectory (some works still needed)

% requires password(s) to be entered if a auth key not available 
% see e.g.:
% https://www.howtogeek.com/66776/how-to-remotely-copy-files-over-ssh-without-entering-your-password/

clear all
close all
% addpath(genpath('/Users/stevenjerjian/Desktop/FetschLab/Analysis/codes/'))
addpath(genpath('C:\Users\yhaile2\Documents\AcademicRelated\CODE_Projects\GitHubCodes\Fetschlab\FLprojects\preliminaries'))



%% decide which files to load

today    = str2double(datestr(now,'yyyymmdd'));

subject = 'lucio';
paradigm = 'dots3DMP'; % grab data for this task paradigm
% paradigms = {'dots3DMPtuning','dots3DMP','RFMapping','VesMapping'};


% dateRange = 20220512:today; % RT 20221215
dateRange = 20220512:20220514; % testing store eyePos and Ideal mp vectors (visTraj-->now written as posTraj & velProfile)
% dateRange = 20230925:today; % RT 20221215
% dateRange = 20200203:20200207;
% dateRange = 20230501:today;
% subject = 'human';
% dateRange = 20190625:20191231; % non-RT
% % dateRange = 20200213:today; % RT

dateStr = num2str(dateRange(1));
for d = 2:length(dateRange)
    dateStr = [dateStr newline num2str(dateRange(d))];
end

% localDir = ['/Users/stevenjerjian/Desktop/FetschLab/PLDAPS_data/' subject '/'];
% localDir = ['/Users/chris/Documents/MATLAB/PLDAPS_data/' subject '/'];
% localDir = ['C:\Users\yhaile2\Documents\CODE\MATLAB\3DMP_Zarya\Zarya_Behavior\3DMP_BehavData\'] ; % YYY %
localDir = ['C:\Users\yhaile2\Documents\AcademicRelated\CODE_Projects\Data\Lucio\Behavior\3DMP\'] ; % YYY % <---- not using subject since folder name is Zarya not 'zarya' as used on server

remoteDir = ['/var/services/homes/fetschlab/data/' subject '/'];

outputPath = ['C:\Users\yhaile2\Documents\AcademicRelated\CODE_Projects\Data\' subject '\Behavior\ConcatBehavSessFiles']; % path to save output concatenated 'data' var to
% to load files directly from mounted homes Volume (be careful not to save
% over to the original file!)
% mountDir  = ['/Volumes/homes/fetschlab/data/' subject '/'];

mountDir = ['Z:\fetschlab\data\' subject '\'] ; % YYY %
% mountDir = ['Z:\fetschlab\data\' subject '\Lucio_NasBehavFiles_Backup\'] ; % YYY % %ollld behav block files

addNexonarDataToStruct = 0; % SJ 08-2021
addDotPositionToStruct = 0; % SJ 01-2022
addMotionTrackingData = 0;
addEyeMovementToStruct = 1;
saveRewardData = 1;         % SJ 06-2022

%% get PDS files from server -- DON'T FORGET VPN
% will skip files that already exist locally, unless overwrite set to 1

useSCP = 0; % 1 - secure copy of files to local folder, 0 - load files directly from mounted drive, save locally only after cleanup
useVPN = 0; % 1 - use proxy VPN (off campus), 0 - use 172 address
overwriteLocalFiles = 1; % set to 1 to always download anew the server copy

getDataFromServer % now also includes pdsCleanup to reduce file size and complexity

% get Nexonar files from server

if addNexonarDataToStruct
    localDir = ['/Users/stevenjerjian/Desktop/FetschLab/PLDAPS_data/' subject '/nexonar/'];
    remoteDir = ['/var/services/homes/fetschlab/data/' subject '/' subject '_nexonar/'] ;

    getDataFromServer % get nexonar data, this will skip pdsCleanup
    
    % re-assign localDirs for createDataStructure
    localDir = ['/Users/stevenjerjian/Desktop/FetschLab/PLDAPS_data/' subject '/'];
    localDirNex = ['/Users/stevenjerjian/Desktop/FetschLab/PLDAPS_data/' subject '/nexonar/'];
end

%% create data structure

createDataStructure

%% optional: save data struct to a mat file so you don't have to repeat the time consuming step

if sum(diff(dateRange)>1)==0
    file = [subject '_' num2str(dateRange(1)) '-' num2str(dateRange(end)) '.mat'];
elseif sum(diff(dateRange)>1)==1
    file = [subject '_' num2str(dateRange(1)) '-' num2str(dateRange(diff(dateRange)>1)) '+' num2str(dateRange(find(diff(dateRange)>1)+1)) '-' num2str(dateRange(end)) '.mat'];
else
    warning('Multiple breaks in date range, could overwrite existing file with different intervening dates!');
    disp('Press a key to continue');
    pause
    file = [subject '_' num2str(dateRange(1)) '---' num2str(dateRange(end)) '.mat'];
end
    
try
data = rmfield(data,'dotPos'); % CAREFUL
catch
end
%%
disp('saving...');
cd(outputPath)
save(file,'data') % YYY %
% save([localDir(1:length(localDir)-length(subject)-1) file], 'data');

% otherwise for larger files will need: 
% save([localDir(1:length(localDir)-length(subject)-1) file], 'data','-v7.3');
% save(file,'data','-v7.3'); % YYY %
% save([localDir file], 'data','-v7.3'); % YYY %

disp('done.');

