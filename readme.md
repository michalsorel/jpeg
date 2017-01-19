## JPEG decompression and denoising with tight frame priors
This repo contains a Matlab code implementing two JPEG decompression algorithms described in the following
two papers. Example images can be found in <http://zoi.utia.cas.cz/jpegrestoration>.

*M. Šorel, M. Bartoš, Efficient JPEG decompression by the alternating direction method of multipliers, ICPR, 2016*  
This conference paper describes the application of the alternating direction method of multipliers (ADMM) 
to the problem of JPEG decompression in the Bayesian MAP formulation with priors given by the 
l1 norm of a tight frame operator. We provide all the necessary math including the proof of an 
expression for the projection on the quantization constraint set. Convergence is extremely fast 
and given the speed of the method, results are competitive with other state-of-the-art methods. 
Note that some recent methods based on neural networks preserve slightly more details 
but can hardly compete in speed.  

*M. Šorel, M. Bartoš, Fast Bayesian JPEG decompression and denoising with tight frame priors, IEEE Transactions on Image Processing, vol. 26, no. 1, pp. 490-501, Jan. 2017*  
This paper provides a more detailed analysis of the MAP formulation of JPEG decompression with 
priors given by the l1 norm of a tight frame operator. Important difference with respect to the 
conference paper is that it considers a noise added in the image before compression. In other words, 
we consider the situation of an image degraded by a noise and then compressed by the JPEG algorithm. 
The task is to recover the original image by maximizing the posterior probability using the same priors 
we use in the JPEG decompression. We also show that using Gaussian approximation instead of 
the hard quantization constraint set surprisingly improves results. We derive how this (different) 
formulation can be efficiently optimized by the ADMM.

All our code is in Matlab. We read and parse JPEG files using the JPEG toolbox of Phil Sallee
that requires compilation to mex files. We provide compiled mex files for 64-bit windows, though.
For regularization, we use primarily the dual-tree complex wavelets (original code is included) but the algorithm 
can use also orthogonal wavelets from the Matlab wavelet toolbox (Haar, db2, db3), learned tight frames of Cai et al. or any combination of them.   

## Files and directories
`start_on_jpeg.m` - elementary example script how to run a reconstruction of a JPEG-file.  
`reconstruct_jpeg.m` - wraper function for reconstruction from JPEG-file. You can use it to 
decompress JPEG files with various parameters. Parameters are described in the file `start_manually.m`.   
`start_manually.m` - script demonstrating various algorithmic options you can set. 

`jpeg_decoder_bregman_color_noise.m` - code of the reconstruction algorithm.  
`dict8all500000.mat`, `dict16all500000.mat` - precomputed learned frames learned on a database of architecture images.  
`proj2QCS.m`, `box_projection.m` - projection functions.  
`additional\` - supporting functions for colorspace transformations and evaluation metric computations (SNR, SSIM).  
`frames\` - functions related to frames  
`images\` - test images  
`jpegtbx_1.4\` - Matlab JPEG Toolbox written by Phil Sallee 9/2003 with mex-files compiled for Win64. For other operating systems see instalation instructions inside.

**Parameters and settings**

The reconstruction algorithm can use different settings. For details see `start_manually.m`, 
section *Reconstruction*. The main options are:

*Method*  
We provide two versions of the JPEG reconstruction algorithm. Rigorous one with
the theoretically correct quantization constraint set described in the conference paper
and one with an approximation of the quantization interval by a Gaussian function
(journal paper).
You can set them as
`par.method_version = 'QCS';` and `par.method_version = 'Gauss';`.

*Image prior, regularization*  
The image prior can be described by the l<sub>1</sub> norm of various tight frames. 
This is specified by the field <em>frame</em> of the parameter structure, see `start_on_jpeg.m`.

* Daubechies wavelets - orthogonal Daubechies wavelets from Matlab Wavelet Toolbox. Requires Wavelet Toolbox. 
Available wavelets are: *Haar, db2, db3*. The fastest option. For example `par.frame={'Haar'};`
* DT-CWT (dual-tree complex wavelets) - better than standard wavelets on diagonal edges. The implementation
we use is probably slower than it could be. Otherwise, it is a good compromise of speed and quality.  
Example: `par.frame={'DT-CWT'};`
* Learned frame - a frame learned on our image database. We use the frame of Cai et al. (see the papers), which
is basically a set of filters (AKA field of experts) enforced to be a Parseval frame. 
You can learn your own frame by dowloading the original code from the web of Jian-Feng Cai
<http://www.math.ust.hk/~jfcai/>.  
Example: `par.frame={'LearnedD'};`

You can also compose larger frames from the basic ones, for example from  *Haar, db2*, and *db3*
by  `par.frame = {'Haar','db2','db3'};`.

*Joint regularization*  
Since the chroma channels in JPEG are subsampled, the most information is stored in the
luminance channel. Our code can use *joint regularization*, an application of the idea of 
group sparsity, which improves details in the chroma channels.  This setting is implicitly on, i.e.
`par.multichannel_regularization = true;`.

*Others*  
There are several other parameters, for example the number of iterations
or the regularization weights. They are documented in the code.


##Terms of Use

This code can be freely used for research purposes by academic organizations.
If you use it in your research, please cite the papers above.