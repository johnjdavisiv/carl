function CA_rules_index = getEnergeticActivity(vm,p2p_window, p2p_threshold, continuity, fs)

%{
Apply the "continuous" and "amplitude" rules for the CARL classifier

John Davis & Marcin Straczkiewicz
14 Oct 2020

Parameters
----------

vm: n x 1 array
        A (possibly filtered) vector magnitude (resultant) acceleration signal, in gravitational
        units (aka g-units). NaN values are permitted but are always returned as 0 (below threshold)

p2p_window: int
        How long of a window, in seconds, should the peak-to-peak amplitude calculations use? CARL
        uses one second windows.

p2p_threshold: float
        What is the minimum peak-to-peak acceleration amplitude (in g) a given window needs in order
        to successfully pass the amplitude rule?

continuity: int
        What is the minimum duration, in seconds, that vm must continuously exceed p2p_threshold in
        order to successfully pass the continuity rule? continuity determines the minimum possible
        duration of running that can be detected. 

fs: int or float
        Sample frequency of vm, in Hz. Same caveats as CARLclassify's fs parameter.


Returns
-------

CA_rules_index: logical array
        A vector the same length as vm that indicates where vm passes the continuity (C) and
        amplitude (A) rules.

Notes
-----

This is an internal function of CARL, so no data validation is performed (all of that happens in the
main call to CARLclassify).

getEnergeticData(my24hr_data,1,1.5,5, 100)
Will return a logical vector that indicates when activity happened that
had a peak-to-peak amplitude of at least 1.5 g in 1-second windows, that
was continuous for at least 5 seconds (in data collected at 100 Hz).    

%}


%% Rule 1: Peak-to-peak amplitude

windowed_r = reshape(vm,round(fs)*p2p_window,[])';

%Calculate peak to peak amplitude for each window
window_p2p = max(windowed_r,[],2) - min(windowed_r,[],2);

p2p_index = repmat(window_p2p > p2p_threshold,[1,round(fs)*p2p_window])';
p2p_index = p2p_index(:);  %Turn into a 1xn logical vector

%% Rule 2: Continous activity
%Where is the logical vector peak_to_peak *not* changing? This is where we
%have 'streaks' of continous values (either 1 or 0)
if nnz(p2p_index) > 0 %If the data has at least one window of energetic activity
    if nnz(p2p_index) == length(vm) %If the input data is ALL above energetic threshold
        streak_index_length = [1, length(vm), 1, length(vm)];
    else
        mydiff = (diff(p2p_index)==0);
        ix=find(mydiff == 0);  % find the indices where these streaks end
        ix(end+1) = length(p2p_index); %Without this you will always miss the LAST streak
        
        if size(ix,2) ~= 1 %only occurs when there is only one streak, just make it a column vector
            ix = ix';
        end

        sl=diff([0; ix]);   %sl = streak length, i.e. how long was this streak?

        for a=1:length(ix) %For each streak that exists...
            ix(a,2) = p2p_index(ix(a,1));
        end

        streak_index_length = [ix(:,1)-sl+1, ix, sl];
    end
    %streak_index_length has four columns
    %The first column tells you the index of the START of a streak
    %The second column tells you the index of the END of a streak. 
    %The third column tells you whether this was a streak of 
    %  below-threshold acceleration (0) or above-threshold acceleration (1). 
    %The fourth column tells you how long (in samples) that streak was.

    % Generate the logical vector to return 
    sxd = zeros(length(vm),1); %sdx = streak day index, index vector to return
    %Apply rule for continuous activity
    for a=1:size(streak_index_length,1) %For each streak...
        if (streak_index_length(a,4) >= continuity*round(fs)) && (streak_index_length(a,3) == 1) 
            %If we had a streak as long as our threshold AND it was a streak of energetic activity... 
            %make a length of ones in rdx that ends at this index
            sxd(streak_index_length(a,2) - streak_index_length(a,4) + 1:streak_index_length(a,2)) = 1;
        end
    end
    CA_rules_index = logical(sxd);
else %If nothing passed the p2p threshold...
    CA_rules_index = false(length(vm),1);
end

end

