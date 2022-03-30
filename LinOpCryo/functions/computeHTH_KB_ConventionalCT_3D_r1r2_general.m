%=========================================================================%
% This function computes the kernel of HTH using the convolution of the
% projection of each basis. It performs the computation for Kaiser-Bessel
% window functions. It computes the convolution using FFT and IFFT, and
% then interpolate it on the grid <k,\theta>.
%-------------------------------------------------------------------------%
%%    Copyright (C)
%     2018 - Biomedical Imaging Group, EPFL
%     Masih Nilchian - masih.nilchian@epfl.ch
%     Laurène Donati - laurene.donati@epfl.ch
%     Citation:  Fast Multiresolution Reconstruction for Cryo-EM, Donati et al.

function y = computeHTH_KB_ConventionalCT_3D_r1r2_general(param)

k1       =   param.k1;
k2       =   param.k2;
k3       =   param.k3;
r1       =   param.r1;
r2       =   param.r2;
scaleRatio=param.scaleRatio;
a        =   param.a;
alpha    =   param.alpha;
m        =   param.m;


% Compute the projection of one basis function

Num_Sample = 1001;
x_value  =   linspace(-2*a,2*a,Num_Sample);
[X,Y]    =   meshgrid(x_value,x_value)    ;
s        =   sqrt(X.^2+Y.^2)                          ; 
X_KBWF   =   zeros(length(x_value),length(x_value));   %X-ray projection of Kaiser_Bessel Window function
z        =   alpha*sqrt(1-(abs(s)/a).^2);
X_KBWF(abs(s)<a) =  (a/besseli(m,alpha))*sqrt(2*pi/alpha)*...
    ((z(abs(s)<a)/alpha).^(m+1/2)).*...
    besseli(m+1/2,z(abs(s)<a));
yConv    =   ifftshift(ifft2(fft2(X_KBWF).*fft2(X_KBWF)))*16*a^2/((Num_Sample-1)^2);


% antialiasing
if (param.aliasingFilter)
    
    yConv = doAntiAliasingFilter2D(yConv, scaleRatio*4*a/(Num_Sample-1)) ;
    
end


% make look-up table
LT       =   yConv((Num_Sample+1)/2,(Num_Sample+1)/2:end);
LTindex  =   x_value((Num_Sample+1)/2:end)  ;

y                =  zeros(2*k1-1,2*k2-1,2*k3-1)                             ;
[K1,K2,K3]       =   meshgrid((-k1+1):(k1-1),(-k2+1):(k2-1),(-k3+1):(k3-1)) ;

for i = 1 : size(r1,2)
    % set the kernel
    if mod(i, 200)==0
        fprintf('computing Hth: i=%d \n', i)
    end
    kInerr1Tmp          =   K1*r1(1,i)+K2*r1(2,i)+K3*r1(3,i)                   ;
    kInerr2Tmp          =   K1*r2(1,i)+K2*r2(2,i)+K3*r2(3,i)                   ;
    KInerAngleTmp       =   sqrt(kInerr1Tmp.^2+kInerr2Tmp.^2)             ;
    y(KInerAngleTmp<(2*a))  =   y(KInerAngleTmp<(2*a))+scaleRatio^4*interp1(LTindex,LT,KInerAngleTmp(KInerAngleTmp<(2*a)));
end
end




