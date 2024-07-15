# deconvTRAP

Code to deconvolve TRAP datasets from their input, with examples.

deconvTRAP identifies the percentage of whole tissue (input) that can be explained by TRAP, i.e. the fraction of the TRAPped cells within the tissue, and the percentage of TRAP that can be explained by the input, i.e. the contamination of the input material that was aspecifically bound by the IP.

Following this, deconvolution is operated and input-free TRAP, as well as TRAP-free input libraries are obtained that can be used for following analyses.

![Illustration of the TRAP deconvolution method and results](https://github.com/ste-depo/deconvTRAP/blob/main/figs/Deconvolution_method.jpg?raw=true)
