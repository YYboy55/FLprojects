function LI = PTE_genLI(datastruct, fieldnames, nback, fieldValues)
% Takes as Input a datastruct, a cell array of fieldnames, vector of nback values, cell array of field values
% Example inputs:
    % datastruct = data;
    % fieldnames = {'modality', 'PDW', 'modality&PDW', 'coherence'};
    % nback = {1, 2, [1, 1], 3};
    % fieldValues = {2, 1, {2, 1}, 0.8};
% 'modality&PDW' here is unique in that two field names are given as a single string.
% This is done when the desired LI output is conditioned on the values in multiple fields.
% Thus, you can see in the third cell of 'fieldValues' 2 values are given to associate with modality and PDW respectively.
% 'nback' allows specification of whether you care about evaluating the data values on this trial (n=0) or past trials (n-1,2,3 etc.)

numConditions = numel(fieldnames);

    % Initialize output variable 'LI'
    if contains(fieldnames{1}, '&')
        subfieldNames = strsplit(fieldnames{1}, '&');
        LI = false(numel(datastruct.(subfieldNames{1})), numConditions);
    else
        LI = false(numel(datastruct.(fieldnames{1})), numConditions);
    end

    % Loop over each fieldname (i) + subfieldName (j)
    for i = 1:numConditions
        fieldName = fieldnames{i};
        fieldValue = fieldValues{i};
        n = nback{i};
        if contains(fieldName, '&')
            % Field name contains '&', split it into individual field names
            subfieldNames = strsplit(fieldName, '&');
            numSubfields = numel(subfieldNames);
            logicalVec = true(size(datastruct.(subfieldNames{1}))); % logicalVec = length(data) rows x 1 column
            
            % Calculate the maximum nback value among subfieldnames within the fieldname
            % maxNback = max(nback{i});

            % Evaluate each subfield against subfieldValue 
            for j = 1:numSubfields
                subfieldName = strtrim(subfieldNames{j});
                subfieldValue = fieldValue(j);
                Nback = n(j);
                subfieldLogicalVec = datastruct.(subfieldName)(((Nback+1):length(datastruct.(subfieldName)))-Nback) == subfieldValue; % True when subfield-specific condition is met

                % Pad the subfieldLogicalVec with false values for Nback rows
                subfieldPadding = false(Nback, 1);
                subfieldLogicalVec = [subfieldPadding; subfieldLogicalVec];
                logicalVec = subfieldLogicalVec;
%                 % Adjust the length of logicalVec if necessary
% Assuming the data from which all subLogiVecs are formed is of equal
% length all subLogiVecs should be same length. Subtract Nback values then
% add back Nback 0s to the vector = take as much as given back = same as
% original length.
%                 subfieldLength = numel(subfieldLogicalVec);
%                 logicalVec = logicalVec(1:subfieldLength) & subfieldLogicalVec;
            end

%             % Pad the logicalVec with false values for maximum nback rows
% This seems unnecessary since padding is occuring for each subfield
% individually, according to it's specific nback value.
%             padding = false(maxNback, 1);
%             logicalVec = [padding; logicalVec];
        else
            % Field name does not contain '&', evaluate it directly
            logicalVec = datastruct.(fieldName)(((n+1):length(datastruct.(fieldName)))-n) == fieldValue; % True when field-specific condition is met

            % Pad the logicalVec with false values for nback rows
            padding = false(n, 1);
            logicalVec = [padding; logicalVec];
        end

        % Store the logical vector in the corresponding column of the output
        LI(:, i) = logicalVec; % LI values < nback+1 remain as 0. These trials are incapable of being assessed for desired past trial effects.
    end
end






