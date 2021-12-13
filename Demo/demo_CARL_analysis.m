%% Showcasing CARL classifier usage on raw free-living accelerometer data

addpath('../CARLclassify function');

% Load up 24hr of 100 Hz freeliving data from the low back: 
load('demo_24hr_actigraph_data_100hz.mat');
%four columns are acceleration (in g-units): x y z r
% (r being resultant, aka "vector magnitude"


%CARL takes vector magnitude input only!
vm = day_xyzr(:,4);
t = linspace(0,24,length(vm));

%Plot daily activity
figure(1);
plot(t, vm, 'color', [0 0 0.5]);
xlim([0 24]);
xlabel('Hours since midnight');
title('Raw data');
ylabel('Acceleration (g)');

%Run CARL, specifying device location, minimum bout length of running (in seconds), and sampling
%frequency (in Hz)
vm_logical_final = CARLclassify(vm, 'torso', 5, 100);
%CARL returns a logical vector the same length as the input vector. 1 = running!


%Set non-running acceleration to NaN, just for plotting purposes
vm_final = vm;
vm_final(~vm_logical_final) = NaN;

figure(2);
hold on;
plot(t, vm, 'color', [0.5 0.5 0.5 0.5]);
plot(t, vm_final, 'color', [0 0 0.5]);
xlim([0 24]);
xlabel('Hours since midnight');
title('Blue: Running. Gray: Non-running')
ylabel('Acceleration (g)');

%Notice how CARL can accurately ingore brief interruptions, as well as the transition to walking at the end
%of the run.
figure(3);
hold on;
plot(t, vm, 'color', [0.5 0.5 0.5 0.5]);
plot(t, vm_final, 'color', [0 0 0.5])
xlim([12.4 12.52]);
title('Close-up on mid-run interruptions');
ylabel('Acceleration (g)');

%% Raw feature extraction


%Show use of CARLfeatures(), which extracts the dominant frequency (in Hz) and peak-to-peak amplitude (in
%g), which are the raw features used by CARL to detect running. 

%Warning: this can take several seconds
[wav_ridge, p2p_val] = CARLfeatures(vm, 100, 16, 8, 1);

figure(4);
scatter(wav_ridge, p2p_val, 20, [0 0 0.5], 'filled', 'markerfacealpha', 0.1);
xlabel('Dominant frequency (Hz)');
ylabel('Peak-to-peak amplitude (g)');
title('Raw features in one-second windows');
