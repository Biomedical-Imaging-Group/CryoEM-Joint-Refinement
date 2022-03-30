%%    Copyright (C) 
%     2018 - Biomedical Imaging Group, EPFL
%     Masih Nilchian - masih.nilchian@epfl.ch	
%     Laurène Donati - laurene.donati@epfl.ch


function R = Euler(phi, theta, psi)
	R = [cos(psi)*cos(theta)*cos(phi)-sin(psi)*sin(phi), cos(psi)*cos(theta)*sin(phi)+sin(psi)*cos(phi), -cos(psi)*sin(theta);
		-sin(psi)*cos(theta)*cos(phi)-cos(psi)*sin(phi), -sin(psi)*cos(theta)*sin(phi)+cos(psi)*cos(phi), sin(psi)*sin(theta);
		sin(theta)*cos(phi), sin(theta)*sin(phi), cos(theta)];
end