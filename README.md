White-Box-LPN
=============

This repository contains all implementations and time samples related to our submission called **LPN-based Attacks in the White-box Setting** to the CHES 2023.

This repository is divided into three directories:

* `implementations` contains the implementations of `timeHDDA.sage`, `timeHODCA.sage`, and `timeWBLPN.sage` used to measure the time of their corresponding algorithms. It also contains our reference implementation of `WBLPN.sage` and `HOWBLPN.sage`. Lastly, this repository contains `processData.sage`, which uses the saved time samples to produce their average and min/max.

* `ourTimeResults` contains all of our time samples used to generate the two figures 1 & 2 of our paper.

* `timeResults` are empty folders that will be filled by time samples by executing `timeHDDA.sage`, `timeHODCA.sage`, and `timeWBLPN.sage` of the `implementations` directory.

All implementations begin with tunable parameters if they exist.

`WBLPN.sage`
------------
This file is the reference implementation of the adaptation of the algorithm "Pooled Gauss" for white-box cryptography. Executing the code will apply WBLPN to random inputs I=20 times, and see if it returns a solution, which in that case would be a false positive. Thereafter, it will apply WBLPN to random inputs that contain a solution for the defined noise rate. If the algorithm does not manage to find one or finds an incorrect one, it prints the number of noised equations of Mpool. Indeed, even if tau =1/4, we cannot be sure that Mpool will not contain more than 1/4 of equations with noise=1.

`HOWBLPN.sage`
--------------
This file is the reference implementation of the higher-order version of WBLPN. Executing the code will run I=20 times HOWBLPN on window size 20 and second order noise rate = 1/2**7, which is the example given on the paper, and print the average. In the paper we did a mistake: we said that it takes 20s, instead of 3s, due to a mistake.

`timeWBLPN.sage`, `timeHODCA.sage`, and `timeHDDA.sage`
------------------------------------------------------
These files contain non-optimized implementations of WBLPN, HODCA, and HDDA; and run them for different order/degree/noise rate and window sizes from 5 to 50, for I=20 iterations. For each of these iterations, the time is measured and recorded in the empty folders within `timeResults`. These implementations resulted in the two figures of our paper.

`processData.sage`
------------------
This file contains a code that opens **our** time measurements, and computes the average, minimum, and maximum for each of them, to print it on the screen. To make this code process your data, simply change every line stating "ourTimeResults" to "timeResults".
