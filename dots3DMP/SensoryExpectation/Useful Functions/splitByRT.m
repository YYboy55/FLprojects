% load lucio_20220512-20230906_clean.mat
% % addpath(genpath('C:\Users\yhaile2\Documents\CODE_Projects\MATLAB\3DMP\SensoryExpectation\Functions_CueWeights_PastTrialEffects\PastTrialEffectsPackage'))
% addpath(genpath('C:\Users\yhaile2\Documents\CODE_Projects\GitHubCodes\Fetschlab\FLprojects\dots3DMP\SensoryExpectation\PastTrialEffectsOnCueW__Package'))

function slowFast_RTidx = splitbyRT(data)
hdng = data.heading;
RT = data.RT;
hdng(hdng>3 & hdng<4) = 3;
hdng(hdng<-3 & hdng>-4) = -3;
hdng(hdng == 1.25) = 1.5;
hdngs = unique(hdng);
m3 = data.modality == 3;

slowRT = false(size(hdng));
fastRT = false(size(hdng));
for h = 1:length(hdngs)

    h_idx = find(hdng==hdngs(h) & m3);
    RT_split = RT(h_idx) > median(RT(h_idx));
    slowTr = h_idx(RT_split); % Calc cue weights for the slower half of mod 3 trials of heading h
    fastTr = h_idx(~RT_split); % Now for the faster half

    slowRT(slowTr) = true; % for a given heading, have marked the slow 50% of trials. GOAL = collect the slow half of trials at each heading and collectively index to these as the slow RT trials = create a dataset of all slow trials :)
    fastRT(fastTr) = true;

end
% Now have LI to fast and slow trials. Can input into cueWcalc 

slowFast_RTidx = [slowRT,fastRT];
end

