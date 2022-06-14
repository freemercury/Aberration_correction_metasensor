# Zernike-Polynomials-MATLAB
Set of functions for 1) quickly generating Zernike polynomials and 2) performing least-squares fits of images using Zernike polynomials.

These functions may be used to quickly generate Zernike polynomials of any radial and azimuthal degree over a circular aperture of any resolution. Additionally, they may be used to perform a quick least-squares fit of any image within a circular aperture using Zernike polynomials, returning the relative coefficients (or "moments", as described by the literature) of each polynomial used in the fit. Many thanks to Chong et. al for proposing a recursive algorithm to calculate the radial portion of the polynomials; it is this algorithm that makes the code as efficient as it is.

Function zernike is used to generate Zernike polynomials. 
Functions zernike_moments and zernike_recreation are used to perform a least-squares fit and recreation of an image using Zernike polynomials.

Please read function descriptions for full instructions on their use.

[![View Quick Zernike polynomial creation and decomposition on File Exchange]
(https://www.mathworks.com/matlabcentral/images/matlab-file-exchange.svg)](https://www.mathworks.com/matlabcentral/fileexchange/72440-quick-zernike-polynomial-creation-and-decomposition)
