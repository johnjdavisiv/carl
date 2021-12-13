function [wavelet_ridge_at_nwin, wavelet_ridge_at_fs1] = getWaveletRidge(vm_raw,fs1,fs2,lp_cutoff,window_size)
%{

Return the dominant frequency of a vector by using the continuous wavelet transform

John Davis & Marcin Straczkiewicz
14 Oct 2020

Parameters
----------

vm_raw: n x 1 array
        A (possibly filtered) vector magnitude (resultant) acceleration signal, in gravitational
        units (aka g-units). Length should be an integer number of seconds (CARL takes care of this
        internally for noninteger-second data) 

fs1: int
        Sample frequency of vm_raw, in Hz. Must be integer (CARL takes care of this but raises warning).


fs2: int
        Sample frequency to downsample vm_filt after filtering. Helps computational speed. Recommend
        using 16 Hz for CARL. More generally, use 2*lp_cutoff.

lp_cutoff: int
        Cutoff frequency for Butterworth lowpass filter. Use -1 to skip filtering. CARL already
        filters so when called internally CARLclassify specifies -1. Generally, use half of fs2.

window_size: int
        Size of windows, in seconds, to average across for calculating dominant frequency (i.e.
        wavelet ridge). Wavelet ridge will be constant across intervals of window_size*fs1. CARL
        uses one second by default.

Returns
-------

wavelet_ridge_at_nwin: float array
        Dominant frequency, in Hz, for the data when split into window_size windows. Use this as
        your feature in activity recogntion. 

wavelet_ridge_at_fs1: float array
        Dominant frequency, in Hz, of the data at its raw sample frequency. Involves no further
        computation than wavelet_ridge_at_nwin (in fact it just uses repelem()). Convenient for
        plotting amplitude data alongside frequency data. 

Notes
-----

The dominant frequency (termed "wavelet ridge") is estimated by averaging the magnitude of the
continuous wavelet transform coefficients across window_size, then finding the maximum. This 
averaging is what permits CARL to detect even very short bouts of running.

Uses analytical Morse wavelet with 4 octaves and 48 voices per octave. Dynamically adjusts if the
input singal is too short, but will raise a warning. This is an internal function of CARL so most
data validation is left to CARLclassify(). Downsampling dramatically increases speed and decreases
memory footprint without any impact on human activity recognition performance. 

%}

%% Lowpass filter + resample + CWT

%Note: if data length (in seconds) is not an integer multiple of the window
%size, the code will skip the remaining data beyond the last full integer
%window.

if min(size(vm_raw)~=1)
    error('Input was not a vector. Did you input a matrix instead?');
end

%Set up constants
nwin = floor(length(vm_raw)/fs1/window_size);
wavelet_ridge_at_nwin = zeros(nwin,1);
norm_cutoff_freq = lp_cutoff/(fs1/2); %set normalized cutoff frequency
filter_order = 2; %Gets doubled to fourth-order in filtfilt()

if lp_cutoff ~= -1
    %Perform lowpass filter
    [bee,ayy] = butter(filter_order,norm_cutoff_freq);  
    vm_filt = filtfilt(bee,ayy, vm_raw);
else %do not filter if lp_cutoff is set to -1
    vm_filt = vm_raw;
end
%     
samp_mean = mean(vm_filt); %Need to detrend to prevent edge effects when resampling
downsamp_vm = resample(vm_filt - samp_mean,fs2,fs1) + samp_mean; %add the sample mean back ('un-detrend')

%Deal with edge cases where signal is very short by dynamically adjusting numOctaves
[minfreq, maxfreq] = cwtfreqbounds(length(downsamp_vm),fs2, 'Wavelet', 'amor');
max_octaves = floor(log2(maxfreq/minfreq));

%Four should be plenty for max
if max_octaves >= 4
    n_octaves = 4;
elseif max_octaves == 3
    n_octaves = 3;
elseif max_octaves == 2
    n_octaves = 2;
else
    n_octaves =1;
    warning('Signal is very short! Only using one octave for CWT. Results may not be accurate.')
end

%Do cwt on downsampled signal
[coeffs, frequencies] = cwt(downsamp_vm, 'amor', fs2,'VoicesPerOctave', 48,'NumOctaves', n_octaves); 
coeffs = abs(coeffs).^2;

%'amor' = analytical morlet wavelet, which should facilitate eventual
%cross-platform use in, e.g. R or Python. Analytical morlet wavelet also
%has the advantage of having an (as the name suggests) analytical central
%frequency.

for i=1:nwin %For each window... (could vectorize if desired)
    this_window_index =((i-1)*fs2+1):((i-1)*fs2+fs2); %Grab this window of data

    %Sum the coefficients across this window
    summed_amplitudes = sum(coeffs(:,this_window_index),2); %sum coefficients across this window
    [pks, locs] = findpeaks(summed_amplitudes); %Get indicies of all the peaks in this window
    [~, max_loc] = max(pks); %find the index of the biggest peak, if it exists (otherwise will be [])

    if ~isempty(locs(max_loc))
        freq_max = frequencies(locs(max_loc));
        %What frequency corresponds to maximum coefficient?
    else
        freq_max = 0; %If NO local maxima (usually when signal flatlines) set to zero
        %Functionally this will almost never be an issue thanks to the continuous
        %amplitude rules. If this is set to NaN instead, you can blow up
        %your logistic classifier! 
    end
    
    wavelet_ridge_at_nwin(i) = freq_max;
end

%Return full vector at raw sample frequency
wavelet_ridge_at_fs1 = repelem(wavelet_ridge_at_nwin,fs1);

end

