function [wav_ridge, p2p_val] = CARLfeatures(vm, fs1,fs2,lp_cutoff,window_size)

%{
% Return the two features used by CARL to identify running. 

John Davis & Marcin Straczkiewicz
14 Oct 2020

Citation
--------

Davis, Straczkiewicz, Harezlak, and Gruber 2020. (Under Review)

Parameters
----------

vm: n x 1 float array
        A column vector of resultant (i.e. Euclidean norm aka vector magnitude) acceleration data
        from a wearable sensor. Row vectors will be automatically converted to a column vector.
        non-vector inputs will raise an error. Units should be in gs (gravitational units).

fs1: int or float
        Sample frequency of original signal, in Hz. 

fs2: int
        Sample frequency in hz for downsampled wavelet transform (use 16 Hz for vanilla CARL classifier)

lp_cutoff: int
        Cutoff frequency for butterworth filter prior to wavelet transform (use 8 Hz for CARL)

window_size: int
        Window size, in seconds. Use 1 second for vanilla CARL settings.

Returns
-------

wav_ridge: float array
        Array of dominant frequency (wavelet ridge) of input data, in Hz. Arranged into
        window_length-second windows.

p2p_val: float array
        Array of peak-to-peak amplitude of filtered input data, in gs. Arranged into
        window_length-second windows.

Notes 
---------

This function is not used by the CARL classifier but is useful for debugging, or for adapting CARL
for your own unique datasets by generating features for a new classifier. 

If vm is not an integer number of seconds, it is reflected at the end to pad the signal.
If fs is not an integer, the algorithm will proceed, but issue a warning and round fs to the
nearest hz. for the purposes of 1s windowing

%}

%% Deal with possible input issues. 
%Note that helper functions assume all of these as given.

%No matrices allowed
if min(size(vm)~=1)
    error('Input was not a vector. Did you input a matrix instead?');
end

%Column vectors only
if size(vm,1) == 1
    vm = vm(:);
end

%Alert to missing values
if sum(isnan(vm)) ~= 0
    warning('%i NaN samples detected in data (%.2f%% of total samples). Interpolating with cubic splines...\n', ...
        sum(isnan(vm)), sum(isnan(vm))/length(vm)*100);
    vm = fillmissing(vm, 'spline');
end

%Warn for fs not integer
if mod(fs1,1) ~= 0
    warning('Noninteger sample frequency of %.2f. Examining data in windows of %i samples. CARL will work just fine but you may want to resample your data to an integer sample frequency prior to calling CARL', ...
        fs1,round(fs1));
end

%expand to integer number of second by reflecting data
%(will truncate at end)

%If the data is not an integer multiple of one second, reflect back points to pad
%(We will truncate at the end using vm_orig_len)
if mod(length(vm),round(fs1)) ~= 0 
    leftover_samples = round(fs1) - mod(length(vm),round(fs1));
    v_end_refl = vm(end-leftover_samples+1:end);
    %Already know vm is column vector
    vm(end+1:end+leftover_samples) = flipud(v_end_refl);
end
   
%% Perform filtering
norm_cutoff_freq = lp_cutoff/(fs1/2); %set normalized cutoff frequency
filter_order = 2; %double second-order filter after zero-lag filtering

%Perform lowpass filter
[bee,ayy] = butter(filter_order,norm_cutoff_freq);  
vm_f = filtfilt(bee,ayy,vm);

%% p2p threshold

windowed_r = reshape(vm_f,round(fs1)*1,[])';


%% Rule 1: Peak-to-peak amplitude

%Calculate peak to peak amplitude for each window
p2p_val = max(windowed_r,[],2) - min(windowed_r,[],2);


%% Generate p2p and cwt ridge features from (filtered) data

%Peak to peak is just max-min in the window

[wav_ridge, ~] = getWaveletRidge(vm_f,round(fs1),fs2, -1,window_size);
%Get the wavelet ridge data: downsample to 16 hz, do not filter (-1), 
%return ridge in 1 sec windows.


end