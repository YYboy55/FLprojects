% GOAL: 
% Arrange trial data in columns, where each column is associated to a single trial heading

% 'heading','behavior' and 'stimulusProperty' must store same size/symmetric vectors or matrices, where each row is a trial.
%  Each column within a variable is treated separately, and considered paired with the same position columns in other variables.

% INPUTS
% 'heading' the heading angle associated to each trial
%'behavior' the behavioral data stored in any indexable field from PDS that is 1s and 0s (or 1s and 2s) e.g 'data.choice'
%'stimulusProperty' is any logical vector, this acts as an index to
%   'behavior', intended use is to index selectively to trials of particular modality, coherence etc.

% OUTPUT:
% 'dataByHeading' is a cell array that stores matrices of trial data sorted
% into columns by heading AND each columns average.
% Rows correspond to the columns of 'behavior'
% Columns: 1) trial vectors of data sorted by heading 2) average of data by heading

% TO DO 
% Ease of use addition: If using same column vector from one variable to pair with multiple
% columns of other vars... e.g. same heading vector for multiple
% behavior dimensions. Add ability to select 'use for all dimensions'
% --> For now just make sure all inputs have equal number of columns, even if that means
% duplicating same column vector.

function dataByHeading = arrangeByHeading(heading,behavior,LogicalIndex) % 'heading','behavior' and 'stimulusProperty' store same size vectors/matrices, where each row is a trial. Each column within a variable is treated separately, and considered paired with the same position columns in other variables.

if nargin < 3
    LogicalIndex = true(size(behavior));
end

if size(heading,2) ~= size(behavior,2)
    heading = repmat(heading,1,size(behavior,2));
end

if size(LogicalIndex,2) ~= size(behavior,2)
    LogicalIndex = repmat(LogicalIndex,1,size(behavior,2));
end

dataByHeading = cell(size(behavior,2),2);
for b = 1:size(behavior,2)
    behav = behavior(:,b);
    stimProp = LogicalIndex(:,b);
    hdng = heading(:,b);
    hdngs = unique(hdng);
    if max(behav) == 2
    behav = behav - 1;
    end
    byHdng = nan(length(hdng(hdng==mode(hdng))),length(hdngs));
    prByHdng = nan(size(hdngs'));
    for h = 1:length(hdngs) % Split along heading
        byHdng(1:length(behav(stimProp & hdng==hdngs(h))),h) = behav(stimProp & hdng==hdngs(h)); % Index to trials of interest and store in single heading column
        prByHdng(1,h) = sum(byHdng(:,h),'omitnan') / sum(~isnan(byHdng(:,h))); % Number of R choices / number of total choices  % % Choice vectors needa be 1 and 0's!
    end
    dataByHeading{b,1} = byHdng; 
    dataByHeading{b,2} = prByHdng;
    % 'dataByHeading' --> Each row is a different behavioral measure (e.g. choice,
    % correct,confidence, etc.) ordered according to order of columns in behavior var.
    % Row 1 of dataByHdng = Column 1 of behavior
    % Columns correspond to 1) byHdng 2) prByHdng
end