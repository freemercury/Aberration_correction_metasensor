# Aberration_correct_Meta_sensor
Title:      Matlab code for "An integrated imaging sensor for aberration-corrected 3D photography"

Version:    1.0 

Copyright:  2022, JIAMIN WU, Yuduo Guo, Chao Deng,  Lu Fang, QIONGHAI DAI

doi:        **"Link" **

Edited based on the reference [1][2][3].

Matlab code for "An integrated imaging sensor for aberration-corrected 3D photography"
==========================================================

This package contains the implementations of the algorithms, including PSF compute, Lightfield data realign, Lens aberration correction, Motion(Turbulence) correction, and Meta-data reconstruction, which are described in detail in the paper: 
JIAMIN WU, YUDUO GUO, CHAO DENG etc, "An integrated imaging sensor for aberration-corrected 3D photography". Nature, 2022.

Please cite our paper if using the code to generate data (e.g., images, tables of processing times, etc.) in an academic publication.

For algorithmic details, please refer to our paper.

Additional Large size Demo data can be downloaded with the following link:
1. Folder: './Data/Calibration/demo.tif':
**"Link"**
   
2. Folder: './Data/Rawdata/ISO_5_5_11.0ms.0.tif':
**"Link"**
----------------
How to use
----------------
The code is tested in MATLAB 2020b (64bit) under the MS Windows 10 64bit version with an Intel i9-10940X CPU, NVIDIA GeForce RTX 2080 Ti GPU and 128GB RAM.

1. Unpack the package
2. Include subdirectory in your Matlab path
3. Run the .m files in './Code/' :
    Meta_PSF_compute.m: Output **Meta PSF**
    Rawdata_realign.m: Input with **Scanning Light field Raw data**
                       To realign LF Raw into Meta images.
    Reconstruction_main.m: Input with **Meta images** & **Meta PSF**
                       To Remove the aberration and reconstruct high resolution images.
    Lens_correction.m: Input with **Calibration Scanning Light field Raw data**
                        To estimate the whole aberration distribution of the Optical lens.
    Dynamic_correction.m: Input with **Meta images** with motion artefacts
                          To eliminate the motion artefacts when shooting the dynamic scence.
    Turbulence_correction.m: Input with **Meta images** with turbulence artefacts
                             To eliminate the motion artefacts when shooting the turbulent scence.
    
 
   a). We have provided some demo data for test of all the main '.m' files described above. Some of the data are too large for GitHub. So we upload them in the Tsinghua cloud drive described before. The link has been realeased above. There are also some origing raw data of the scanning light field which are included in the Link below. These are the relevant data of the published paper. All of them can be processed through the pipeline with some hyper parameters fine tuned.
   b). The code for 3D photography of Meta images are refereced from paper[3]. The origin code for test can be downloaded here "https://github.com/hotndy/LFDepth_POBR.git"

* Download more Raw scanning light field data from the following link.
https://cloud.tsinghua.edu.cn/d/bb4610b6e82349e89caa/

----------------
Citation 
---------------- 
If you use this code please cite the companion paper where the orginal method appeared:

Wu, J., Guo, Y., Deng, C. et al. "An integrated imaging sensor for aberration-corrected 3D photography". Nature (2022). **"DoiLink"**


----------------
IMPORTANT NOTE 
---------------- 
Should you have any questions regarding this code and the corresponding results, please contact Yuduo Guo (gyd20@mails.tsinghua.edu.cn).

Reference:
1.  R. Prevedel, Y.-G. Yoon, M. Hoffmann, N. Pak, G. Wetzstein, S. Kato, T. Schr?del, R. Raskar, M. Zimmer, E. S. Boyden, and A. Vaziri, 
     "Simultaneous whole-animal 3D imaging of neuronal activity using light-field microscopy," Nature Methods, 2014, 11(7): 727-730.
2.  Lu Z, Wu J, Qiao H, et al. "Phase-space deconvolution for light field microscopy," Optics express, 2019, 27(13): 18131-18145.
3.  J. Chen, J. Hou, Y. Ni, and  L.-P. Chau, Accurate light field depth estimation with superpixel regularization over partially occluded regions. IEEE Transactions on Image Processing.

