%  g = KaiserBesselProjection(m, alpha, a, s)
% 
%  X-ray transform of the kaiser bessel window function
% 
%  Input arguments:
%  m    :    smoothness parameter of KBWF.
%  alpha:    window taper parameter.
%  a    :    support of the KBWF.
%  s    :    s should be positive e.g. 0:.1:a.
%  ---------------------------------------------------------
%  Output arguments:
%  g  :    1D function. X-ray transform of the KBWF along one radius. Since it is isotropic function, this 1D information is sufficient.
%  ---------------------------------------------------------
%  Author: Masih Nilchian
%  Date  : 1 OCT 2015
%  ---------------------------------------------------------
%  example:
%  g   = KaiserBesselProjection(2,8,2,0:.1:2);
%  ---------------------------------------------------------

%%    Copyright (C) 
%     2018 - Biomedical Imaging Group, EPFL
%     Masih Nilchian - masih.nilchian@epfl.ch	
%     Laurène Donati - laurene.donati@epfl.ch

function p = KaiserBesselProjection(m, alpha, a, s)
	tmp = sqrt(1 - (s/a).^2);
	p = a ./ besseli(m, alpha) .* sqrt(2*pi/alpha) .* tmp.^(m+0.5) .* besseli(m+0.5, alpha*tmp);
