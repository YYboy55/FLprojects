function [numVOx2_byCohXdelta,numVOx2_byDelta,avgPctToI,pctTrVO2,vesOx2perBlck,cohs] = numVOx2(data) % numVOx2(filename,file_location)

addpath('C:\Users\yhaile2\Documents\CODE_Projects\MATLAB\3DMP\3DMP_functions')
addpath('C:\Users\yhaile2\Documents\CODE_Projects\MATLAB\3DMP\SensoryExpectation\Functions_CueWeights_PastTrialEffects\PastTrialEffectsPackage')

% if nargin < 2
%     file_location = 'C:\Users\yhaile2\Documents\CODE_Projects\3DMP_Data\BehaviorDataFiles_L\CleanDataFiles_L'; % Default file location if not specified as input
% end
% addpath(file_location)
% %%
% load(filename,'data')

uniq_fnames = unique(data.filename);
TF_VesOx2 = PTE_genLI(data,{'modality&modality&modality'},{[0,1,2]},{[3,1,1]}); % data, fieldnames, nback, fieldvalues 
for b = 1:length(uniq_fnames) % for each block
    [b_length, b_locat] = measureBlckLength(uniq_fnames{b},data.filename);
    vesOx2perBlck(b) = sum(b_locat & TF_VesOx2); % Num of trials in this block that are vesOx2
    b_lengths(b) = b_length;
    pctTrVO2(b) = vesOx2perBlck(b) / b_lengths(b);
end

avgPctToI = mean(pctTrVO2(b_lengths > 500)); % ~ 1% are Trials of Interest before including enough extraVes reps to match num comb reps
% numVOx2 = sum(vesOx2perBlck);


if length(unique(data.delta))>3
    data.delta(data.delta<0) = -3;
    data.delta(data.delta>0) = 3;
end
deltas = unique(data.delta);
cohs = unique(data.coherence);


for d = 1:length(deltas)
    numVOx2_byDelta(1,d) = sum(data.delta==deltas(d) & TF_VesOx2); % Total number of ToI for each delta
    for c = 1:length(cohs)
        numVOx2_byCohXdelta(c,d) = sum(data.delta==deltas(d) & data.coherence==cohs(c) & TF_VesOx2); % Total number of ToI for each delta
    end
end

end


% TF_VesOx2 = PTE_genLI(data,{'modality&modality&modality'},{[0,1,2]},{[3,1,1]}); % data, fieldnames, nback, fieldvalues
% [vesW, vesW_avg] = PTE_cueWcalc(data,TF_VesOx2,300,10)


