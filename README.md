# Data generator with a control over linear separability

We propose a new data generator that is useful for a systematic benchmarking of algorithms for classification and clustering.

## Features
* A user can adjust:
  * how many pairs of classes must be linearly non-separable
  * the number of classes
  * the number of data-points inside a class
  * the probability distribution of data-points
  * the minimal distance between each pair of classes
  * the shape of a point-set that forms a class
* 38 different shapes of classes of various difficulty levels are available.
* The output is a two-dimensional dataset.
* It is easy to use generator in a batch mode by calling the function `createDataset()` with different parameters.

## Getting started
See `examples.m` for some demonstrational examples.

## Publications
* [Nejc Ilc, "Clustering Based on Weighted Ensemble," PhD thesis, University of Ljubljana, 2016](http://dx.doi.org/10.13140/RG.2.2.12151.62882).
* (to appear in) Nejc Ilc, "Data generator with a control over linear separability," Proceedings of 26th International Electrotechnical and Computer Science Conference, Portorož, Slovenia, September 2017. Ljubljana: IEEE Region 8, the Slovenian section of IEEE.

## Acknowledgements
In this project we reused the code from:
* Michael Chen (`sqdistance2.m`)
* Nicolo Giorgetti and Niels Klitgord ([A Matlab MEX Interface for the GLPK library](http://glpkmex.sourceforge.net/))
