% Specify what modalities are present in a given block and index to blocks
% that meet a particular criteria

% Index to blocks by their unique filenames
% Iterate through each block and check what modalities were used for that block.

function [TrIdx,desired_blockIndex,block_names] = blockTypeIndex(data,mod_selection) % ,coh_selection

% Identify which blocks contain and exclude particular trial features (modality, coherence etc.)
% First column of 'mod_selection' is desired mods, column 2 is mods to exclude from block

% Outputs:
% An index to all trials belonging to blocks that meet selection criteria (TrIdx)
% An index to the names of blocks that meet criteria (desired_blockIndex)
% And the corresponding list of all block names to which block index applies (block_names)


block_names = unique(data.filename);
fnames_trList = data.filename;
desired_blockIndex = false(size(block_names));
TrIdx = false(size(data.filename));

for u = 1:length(block_names)
    currBlock = block_names(u);
    blck_locat = strcmp(fnames_trList,currBlock); % make index to trials labeled with 'uniq_fname'
    block_mods = unique(data.modality(blck_locat)) ;   % Now ask what mods are present for this block
    if size(mod_selection,2) == 1
        desired_mods = mod_selection(:,1);
        if length(desired_mods) == length(block_mods) && all(block_mods == desired_mods)
            desired_blockIndex(u) = true; 
        end
    else
        desired_mod = mod_selection(:,1);
        excluded_mod = mod_selection(:,2);
%         if length(desired_mod) > 1 && all(block_mods == desired_mods) &&
%         all(all(block_mods~=excluded_mod))  % Make so treat different if
%         multiple desired-mods are given, in addition to a excluded_mod

        if any(block_mods==desired_mod) && all(all(block_mods~=excluded_mod)) % Are curr block mods compatible with desired and excluded mods?
            desired_blockIndex(u) = true;
        end
    end
    if desired_blockIndex(u) == true
        TrIdx(blck_locat) = true;
    end
end

end


% Selecting blocks by more than just mod composition.. e.g coh present in block
% if nargin > 1
%     desired_cohs = coh_selection(1,:);
%     excluded_cohs = coh_selection(2,:);
% 
%     
%     block_cohs = unique(data.coherence(blck_locat)) ;
