# CARL: a running recognition algorithm for free-living accelerometer data

The CARL classifier is an algorithm built to identify and extract bouts of running from raw, free-living accelerometer data. It can be used for acceleration data collected at the wrist or anywhere on the torso.

The CARL classifier is currently implemented in MATLAB. You will need the signal processing and wavelet toolboxes. Python and R implementations are coming "soon:tm:."

You can find the data accompanying our paper here: **FIGSHARE_DOI**

If you use our data, or the CARL classifier itself, please [cite our manuscript](https://iopscience.iop.org/article/10.1088/1361-6579/ac41b8):  
* Davis, Straczkiewicz, Harezlak, and Gruber. CARL: A running recognition algorithm for free-living accelerometer data. *Physiological Measurement*.	**doi:** 10.1088/1361-6579/ac41b8

## The CARLclassify() function

**CARLclassify.m**  Detect bouts of running in vector-magnitude (resultant) acceleration data from a wearable sensor.  

### Parameters
----------

* **vm**: *n x 1 float array*

A column vector of resultant (i.e. Euclidean norm aka vector magnitude) acceleration data from a wearable sensor. Row vectors will be automatically converted to a column vector. Non-vector inputs will raise an error. Units should be in gs (gravitational units).

* **device_location**: *string*

Either 'torso' or 'wrist'. Use 'torso' for data from the low back, hip, waist, chest, or head. Use 'wrist' for data from the wrist or forearm. 'torso' will probably work best for data from the upper arm, but you may want to try both. CARL does not currently support data collected at the foot ankle, shin, or thigh.

* **continuity**: *int*

Minimum run bout duration, in seconds. CARL will ignore any bouts of running that areshorter than (continuity) seconds. Bouts as short as three seconds are supported. CARL will continue even with shorter bouts, but will raise a warning because of the loss of precision in estimating the dominant frequency. Denoted $\tau$ in our paper.

* **fs**: *int or float*

Sampling frequency of original signal, in Hz. CARL has been validated on data ranging from 20 to 256 Hz and should work for arbitrarily high sampling frequencies. CARL will raise a warning if fs is not an integer, and the analysis will proceed using data windowed to the nearest integer frequency. In most cases noninteger sampling frequencies are not a problem.


### Returns
-------

* **vm_final_logical**: *logical*

A vector of 0s and 1s. Will be equal to 0 at timepoints where no running was detected; will be equal to 1 during bouts of running. Note that the 1s will come in bouts no shorter than fs*continuity.


### Notes
---------

CARL has been validated on both free-living and in-lab data from a wide variety of devices, subjects, and device locations. It uses the resultant acceleration amplitude and dominant frequency, computed using the continuous wavelet transform, to identify bouts of continuous running. In most scenarios, CARL shows excellent accuracy.

Activities that tend to cause false positives are typically those that generate amplitude and frequency content similar to that of running. In our validation work, jump-roping and elliptical machine use were the two primary sources of false positives. In rare cases rapidly ascending or descending stairs very quickly may also be flagged as running.

False negatives may occur with very slow running, or (when using wrist-worn sensors) when a subject grabs on to a siderail for support when running on a treadmill. Our validation work also noted that data collected on loosely-secured smartphones tended to be identified with less accuracy, likely because of the large mass of the sensor and a poor physical coupling between the sensor and the body

Please cite our work if you use the CARL classifier or its supporting data in your research. If you encounter any issues or bugs, please submit a GitHub issue or contact the authors.


## Basic usage

`vm_logical_final = CARLclassify(vm, 'torso', 5, 100);`
The above line runs the CARL classifier, indicating our accelerometer was on on the torso


*See also the `demo_CARL_analysis.m` script in the `\Demo` folder*

The CARL classifier operates on vector-magnitude acceleration (aka resultant acceleration or Euclidean norm), expressed **in *g*-units**. Calculating vector magnitude in MATLAB is easy:

`vm = sqrt(x.^2 + y.^2 + z.^2)  # assuming x, y, z are vectors`

You will need the files in the `/CARLclassify function/` directory on your working path.


## Caveats

CARL will NOT work out of the box with AC-response accelerometers, of the type sometimes used to measure impact shock in running injury research. These sensors cannot detect very low frequency accelerations; the reaction force from the ground (at +1 g) is one such "low frequency acceleration." In fact, its frequency is zero Hz!  

Virtually all wearable accelerometers are DC-response accelerometers, which work just fine with CARL. As a hack fix, if you have AC-response data and you know which way the device was oriented (and you know it stayed mostly oriented the same direction), just add +1 g to the positive vertical axis of the accelerometer to "fake" the presence of gravity. Then just use the logical vector to index into the original data.

## Other included functions

List
* of all the functions
* Included in the folder
* and how you might use them
