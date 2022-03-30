
% Compute the 3D KBWF

% =========================================================================
%
% Copyright 2010-2017 Biomedical Imaging Group at the EPFL.
% Author: masih.nilchian@epfl.ch
% Conditions of use: You are free to use this code for research or
% educational purposes. In addition, we expect you to include adequate
% citations and acknowledgments whenever you present or publish results
% that are based on it.
% Reference: Donati, Laurene, et al. "Compressed sensing for STEM
%            tomography." Ultramicroscopy 179 (2017): 47-56.
%
% =========================================================================

function  Y  =  KaiserBesselWindowFunction3D(x,alpha,a,m)

[X1,X2,X3]  =   ndgrid(x);
X           =   sqrt(X1.^2+X2.^2+X3.^2)    ;

Y           =   zeros(size(X));
z           =   alpha*sqrt(1-(abs(X)/a).^2);

Y(abs(X)<a) =  (1/besseli(m,alpha))*((z(abs(X)<a)/alpha).^m).*...
                           besseli(m,z(abs(X)<a));
