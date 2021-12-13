function [n_bouts, bout_start_ind, bout_end_ind, bout_length] = getBouts(v_logical)
%{
Count "bouts" (uninterrupted streaks) of nonzero values in a logical vector

John Davis & Marcin Straczkiewicz
14 Oct 2020

Parameters
----------

v_logical: n x 1 logical array
        Array of 0s and 1s (false and true) in which you are seeking bouts of 1 (true) values.

Returns
-------

n_bouts: int
        Number of bouts of true (1) values in the array. Will be zero if no bouts exist.

bout_start_ind: n_bouts x 1 array
        The ith element contains the start index of the ith bout of nonzero values in v_logical.
        Will be empty [] if n_bouts == 0.

bout_end_ind: n_bouts x 1 array
        The ith element contains the end index of the ith bout of nonzero values in v_logical. Will
        be empty [] if n_bouts == 0;

bout_length: n_bouts x 1 array
        The ith element contains the length, in samples, of the ith bout. Trivial but useful. Will
        be empty [] if n_bouts == 0

Notes 
---------

This is a helper function that is extremely useful for extracting the raw data from bouts of running
detected using the CARL classifier. Using a for-loop, you can loop through all n_bouts of running
and isolate the data using the start and end indices. See the demo on the CARL GitHub page.

%}

if min(size(v_logical)) ~= 1
    error('getBouts only accepts vector inputs.');
end

%Transpose vector to avoid problems
if size(v_logical,2) ~= 1
    v_logical = v_logical'; %make sure it is a column vector
end


if nnz(v_logical) == 0 %If NO activity at all...just return all NaNs  
    n_bouts = 0;
    bout_start_ind = [];
    bout_end_ind = [];
    bout_length = [];

elseif nnz(v_logical) == length(v_logical) %If the input data is ONE streak of all 1s
    n_bouts = 1;
    bout_start_ind = 1;
    bout_end_ind = length(v_logical);
    bout_length = length(v_logical);

else 
    %Find streaks of continuous 1s
    mydiff = (diff(v_logical)==0);
    ix=find(mydiff == 0);  % find the indices where these streaks end
    ix(end+1) = length(v_logical); %Without this you will always miss the LAST streak
    if size(ix,2) ~= 1 %only occurs when there is only one streak
        ix = ix'; %transpose for consistency
    end      

    sl=diff([0; ix]);     % sl = streak length, i.e. how long was this streak?

    for a=1:length(ix) %This is not vectorized but not a big deal
        ix(a,2) = v_logical(ix(a,1));
    end

    streak_index_length = [ix(:,1)-sl+1, ix, sl]; 

    %streak_index_length (and one_streaks) have four columns:
    %The first and second tell you the index of the START and END of a streak. 
    %The third column tells you whether that was a streak of EC content (1) or non-EC (0)
    %The fourth column tells you how long (in samples) that streak was.

    one_streaks = streak_index_length(streak_index_length(:,3) == 1,:);

    n_bouts = size(one_streaks,1);
    bout_start_ind = one_streaks(:,1);
    bout_end_ind = one_streaks(:,2);
    bout_length = one_streaks(:,4);
    
end

end

