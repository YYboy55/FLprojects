function [hdg,modality,coh,delta,ntrials] = dots3DMP_create_trial_list(hdgs,mods,cohs,deltas,nreps)

if nargin < 5, nreps  = 200; end
if nargin < 4, deltas = [-3 0 3]; end
if nargin < 3, cohs   = [0.1 0.5 0.9]; end
if nargin < 2, mods   = [1 2 3]; end
if nargin < 1, hdgs   = [-10 -5 -2.5 -1.25 -eps eps 1.25 2.5 5 10]; end

%% build trial list
% (can't just randsample the above vectors, because certain combinations of
% modality, coh, delta etc are invalid)
    
% repeat heading list once for ves, ncoh for vis, and ncoh x ndelta for comb
numHdgGroups = any(ismember(mods,1)) + ...
               any(ismember(mods,2)) * length(cohs) + ...
               any(ismember(mods,3)) * length(cohs)*length(deltas);
hdg = repmat(hdgs', numHdgGroups, 1);

lh = length(hdgs);
coh = nan(size(hdg));
delta = nan(size(hdg));
modality = nan(size(hdg));

% kluge for ves, call it the lowest coh by convention
if any(ismember(mods,1))
    coh(1:lh) = cohs(1); 
    delta(1:lh) = 0;
    modality(1:lh) = 1;
    last = lh;
else
    last = 0;    
end

if any(ismember(mods,2)) % loop over coh for vis
    for c = 1:length(cohs)
        these = last+(c-1)*lh+1 : last+(c-1)*lh+lh;
        coh(these) = cohs(c);
        delta(these) = 0;
        modality(these) = 2;
    end
    last = these(end);
end

if any(ismember(mods,3)) % loop over coh and delta for comb
    for c = 1:length(cohs)
        for d = 1:length(deltas)
            here = last + 1 + (c-1)*lh*length(deltas) + (d-1)*lh;
            these = here:here+lh-1;
            coh(these) = cohs(c);
            delta(these) = deltas(d);
            modality(these) = 3;
        end
    end
end

% now replicate times nreps and shuffle (or not):
condlist = [hdg modality coh delta];
% trialTable = repmat(condlist,nreps,1); 
trialTable = Shuffle(repmat(condlist,nreps,1),2);
    % why shuffle? well, sometimes it's easier to find particular trial
    % types to test out when debugging
hdg      = trialTable(:,1);  
modality = trialTable(:,2);  
coh      = trialTable(:,3);  
delta    = trialTable(:,4);  
ntrials  = length(trialTable);



