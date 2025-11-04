function [blck_length,blck_locat] = measureBlckLength(uniq_fname,filenames_trlist)
% determine number of trials in a specific block or blocks. Finds number of times a particular filename is used.
% uniq_fname is a single filename for the block of interest
% filenames is the listing of all trials filenames, aka. 'data.filename'
blck_locat = strcmp(filenames_trlist,uniq_fname); % finds trials labeled with 'uniq_fname'
blck_length = sum(blck_locat);
end

