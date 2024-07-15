# deconvTRAP

Code to deconvolve TRAP datasets from their input, with examples.

deconvTRAP identifies the percentage of whole tissue (input) that can be explained by TRAP ($\alpha$), i.e. the fraction of the TRAPped cells within the tissue, and the percentage of TRAP that can be explained by the input ($\beta$), i.e. the contamination of the input material that was aspecifically bound by the IP.

Starting from TRAP and input libraries (Fig A), deconvolution is operated: input-free TRAP, as well as TRAP-free input libraries are obtainedi (Fig A').

![fig1](figs/Deconvolution_method.jpg?raw=true)

An automatic procedure identifies contamination fractions $\alpha$ and $\beta$ (Fig B) and was able to recapitulate endothelial abunance within different mouse tissues (Fig C). When deconvTRAP was applied to a dataset of endothelial cells (Fig D), it was able to separate endothelial cells from their input in the PCA space (Fig D').

![fig2](figs/Deconvolution_results.jpg?raw=true)
