function vm_final_logical = CARLclassify(vm,device_location, continuity, fs)

%{
Detect bouts of running in vector-magnitude (resultant) acceleration data from a wearable sensor.

CARL: Continuous Amplitude Running Logistic classifier
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
        Non-vector inputs will raise an error. Units should be in gs (gravitational units).

device_location: string
        Either 'torso' or 'wrist'. Use 'torso' for data from the low back, hip, waist, chest, or
        head. Use 'wrist' for data from the wrist or forearm. 'torso' will probably work best for
        data from the upper arm, but you may want to try both. CARL does not currently support data
        collected at the foot ankle, shin, or thigh.

continuity: int
        Minimum run bout duration, in seconds. CARL will ignore any bouts of running that are
        shorter than (continuity) seconds. Bouts as short as three seconds are supported. CARL will
        continue even with bouts as short as two or one second, but will raise a warning because of
        the loss of precision in estimating the dominant frequency. Denoted $\tau$ in our paper.

fs: int or float
        Sampling frequency of original signal, in Hz. CARL has been validated on data ranging from 20
        to 205 Hz and should work for arbitrarily high sample frequencies. CARL will raise a warning
        if fs is not an integer, but the analysis will proceed using data windowed to the nearest
        integer frequency. In most cases noninteger sample frequencies are not a problem. 


Returns
-------

vm_final_logical: logical
        A vector of 0s and 1s. Will be equal to 0 at timepoints where no running was detected; will
        be equal to 1 during bouts of running. Note that the "1"s will come in bouts no shorter than
        fs*continuity.


Notes 
---------

CARL has been validated on both free-living and in-lab data from a wide variety of devices,
subjects, and device locations. It uses the resultant acceleration amplitude and dominant frequency,
computed using the continuous wavelet transform, to identify bouts of continuous running. In most
scenarios, CARL shows excellent accuracy. 

Activities that tend to cause false positives are typically those that generate amplitude and 
frequency content similar to that of running. In our validation work, jump-roping and elliptical 
machine use were the two primary sources of false positives. In rare cases rapidly ascending or 
descending stairs very quickly may also be flagged as running. 

False negatives may occur with very slow running, or (when using wrist-worn sensors) when a subject
grabs on to a siderail for support when running on a treadmill. Our validation work also noted that
data collected on loosely-secured smartphones tended to be identified with less accuracy, likely
because of the large mass of the sensor and a poor physical coupling between the sensor and the body

Please cite our work if you use the CARL classifier or its supporting data in your research. If you
encounter any issues or bugs, please submit a GitHub issue or contact the authors

%}

%Load logistic parameters
load('CARL_final_logistic_parameters_02_03_2021.mat', 'back_final_theta', ...
    'back_best_p', 'wrist_final_theta', 'wrist_best_p');

%Set parameters based on location
switch lower(device_location) 
    case 'torso'
        p2p_threshold = 1.0; %g, comes from HARE2 data
        theta = back_final_theta;
        threshold = back_best_p;
    case 'wrist'
        p2p_threshold = 1.0; %g, comes from HARE2 data
        theta = wrist_final_theta;
        threshold = wrist_best_p;
    otherwise
        error('Invalid location. Specify ''wrist'' or ''torso''');
end


%% Deal with possible input issues. 
%Note that helper functions assume all of these as given.

%Check for possibly erroneous signals (wrong units, wrong accelerometer type)
if mean(vm) > 7.0
    warning('Acceleration data had a mean of %.2f. Are you sure that your data are in g-units?', mean(vm));
end

%No matrices allowed
if min(size(vm)~=1)
    error('Input was not a vector. Did you input a matrix instead?');
end

%Column vectors only
if size(vm,1) == 1
    vm = vm(:);
end

%Input too short
if length(vm) < round(fs)*continuity
    error('vm is length %i but your desired continuity length expects bouts of at least %i samples of running',...
        length(vm), round(fs)*continuity);
end

%Alert to missing values - Some wireless devices drop a sample here or there. 
if sum(isnan(vm)) ~= 0
    warning('%i NaN samples detected in data (%.2f%% of total samples). Interpolating with cubic splines...\n', ...
        sum(isnan(vm)), sum(isnan(vm))/length(vm)*100);
    vm = fillmissing(vm, 'spline');
end

%Warn for fs not integer.
if mod(fs,1) ~= 0
    warning('Noninteger sample frequency of %.2f. Examining data in windows of %i samples. CARL will work just fine but you may want to resample your data to an integer sample frequency prior to calling CARL', ...
        fs,round(fs));
end

%expand to integer number of seconds by reflecting data (will truncate again at end)
vm_orig_len = length(vm);

%If the data is not an integer multiple of one second, reflect back points to pad at the end
%(We will truncate at the end using vm_orig_len)
if mod(length(vm),round(fs)) ~= 0 
    leftover_samples = round(fs) - mod(length(vm),round(fs));
    v_end_refl = vm(end-leftover_samples+1:end);
    %We already know vm is column vector
    vm(end+1:end+leftover_samples) = flipud(v_end_refl);
end
   
%% Perform lowpass filtering

downsamp_freq = 16; %Hz
lp_cutoff = 8; %Lowpass filter cutoff - Hz
norm_cutoff_freq = lp_cutoff/(fs/2); %set normalized cutoff frequency
filter_order = 2; %double second-order filter after zero-lag filtering

[bee,ayy] = butter(filter_order,norm_cutoff_freq);  
vm_f = filtfilt(bee,ayy,vm);

%% Apply continuous amplitude criteria 
%set to NaN places that do not have amplitudes
%above ec_p2p for at least ec_continuity seconds

vm_EC = vm_f; %Don't disturb original filtered data
vm_EC_logical = getEnergeticActivity(vm_f,1,p2p_threshold,continuity,fs); %Store results of getEnergeticActivity as logical vector
vm_EC(~vm_EC_logical) = NaN; %Set 1sec windows that did not pass EC rules to NaN

%% "Squash" the energetic activity into one continuous vector to compute features
%Note that this vector will always be an integer number of seconds long

if nnz(vm_EC_logical) == 0 %If NO energetic activity at all...just return all 0
    vm_final_logical = false(length(vm),1);   
    return %Exit the function immediately - MUCH faster!
elseif nnz(vm_EC_logical) == length(vm) %If the input data is ALL above energetic threshold, all ones
    streak_index_length = [1, length(vm), 1, length(vm)];
else  %This will be the usual case - some possible running, some not
    %Find streaks of continuous activity
    mydiff = (diff(vm_EC_logical)==0);
    ix=find(mydiff == 0);  % find the indices where these streaks end
    ix(end+1) = length(vm); %Without this you will always miss the LAST streak
    
    if size(ix,2) ~= 1 %only occurs when there is only one streak
        ix = ix';
    end      

    sl=diff([0; ix]);
    %streak length, i.e. how long was this streak?

    for a=1:length(ix)
        %This is not vectorized but it is not a big deal unless you have a
        %crazy huge number of bouts (in which case you should chunk out your signal, obviously)
        ix(a,2) = vm_EC_logical(ix(a,1));
    end

    streak_index_length = [ix(:,1)-sl+1, ix, sl]; 
end

EC_streaks = streak_index_length(streak_index_length(:,3) == 1,:);
%streak_index_length (and EC_streaks) have four columns:
%The first and second tell you the index of the START and END of a streak. 
%The third column tells you whether that was a streak of EC content (1) or non-EC (0)
%The fourth column tells you how long (in samples) that streak was.

%Set up data for indexing and squashing
n_EC_samples = sum(EC_streaks(:,4)); %How many datapoints passed the EC rules?
EC_seg_end_ind = cumsum(EC_streaks(:,4));
EC_seg_start_ind = EC_seg_end_ind - EC_streaks(:,4) + 1;
squash_EC_data = nan(n_EC_samples, 1); %squash the EC data together 

%Again not vectorized but this is far from the bottleneck in this code
for a=1:size(EC_streaks,1) %For each bout of EC-passing activity...
    %Put it into the squashed window
    data_ind = EC_streaks(a,1):EC_streaks(a,2);
    squash_ind = EC_seg_start_ind(a):EC_seg_end_ind(a);
    squash_EC_data(squash_ind) = vm_f(data_ind);
end

%% Generate p2p and cwt ridge features from (filtered) data

squash_window = reshape(squash_EC_data,round(fs),[])'; %Window data, 1sec windows
squash_p2p = max(squash_window,[],2) - min(squash_window,[],2); 
%Peak to peak is just max-min in the window

[squash_cwt_ridge, ~] = getWaveletRidge(squash_EC_data,round(fs),downsamp_freq, -1,1);
%Get the wavelet ridge data: downsample to 16 hz, do not filter (-1), 
%return ridge in 1 sec windows.

%% Apply the logistic classifier to predict what is running
%(only within the squashed data that passed the EC rules the first time)

%Run RLR depending on location
switch lower(device_location)  
    case 'torso'
        X = [ones(size(squash_cwt_ridge,1),1), squash_cwt_ridge, squash_p2p];
        LR_predict_windows = predictLogistic(theta, X, threshold); %See bottom for this fxn
        LR_predict_squash_vector = repelem(LR_predict_windows,round(fs));
        
    case 'wrist' %Remember wrist has a quadratic feature (CWT ridge^2)
        X = [ones(size(squash_cwt_ridge,1),1), squash_cwt_ridge, squash_cwt_ridge.^2, squash_p2p];
        LR_predict_windows = predictLogistic(theta, X, threshold); %See bottom for this fxn
        LR_predict_squash_vector = repelem(LR_predict_windows,round(fs));
    otherwise
        error('Invalid device location!');
end

%Now put back the squashed data to its original place
vm_LR_logical = zeros(length(vm_f),1);

for a=1:size(EC_streaks,1) %For all bouts of running...
    data_ind = EC_streaks(a,1):EC_streaks(a,2);
    squash_ind = EC_seg_start_ind(a):EC_seg_end_ind(a);
    
    %take LR logical from squash index and put it in into the RLR
    %indicator at the data_index location
    vm_LR_logical(data_ind) = LR_predict_squash_vector(squash_ind);
end

%% Set periods of nonrunning to NaN

%Remove nonrunning from vm_EC
vm_EC_LR = vm_EC; %save this part
vm_EC_LR(~logical(vm_LR_logical)) = NaN;


%% Apply continuity rules again and return final values

vm_EC_LR_cont_logical = getEnergeticActivity(vm_EC_LR,1,p2p_threshold,continuity,fs);
vm_final_logical = logical((vm_EC_logical & vm_LR_logical & vm_EC_LR_cont_logical));

%Truncate if we had to pad
vm_final_logical = vm_final_logical(1:vm_orig_len);

end

%Few quick utility functions for logistic reg
function p = predictLogistic(theta, X, threshold)
pred = sigmoid(X*theta);
p = (pred > threshold);
end

function g = sigmoid(z)
g = 1.0 ./ (1.0 + exp(-z));
end

