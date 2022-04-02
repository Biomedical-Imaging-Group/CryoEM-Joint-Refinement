classdef LinOpPBTShift <  LinOp
    % TODO: Put here the description of your LinOp
    %
    % :param parName: DESCRIPTION
    %
    % All attributes of parent class :class:`Cost` are inherited.
    %
    % **Note**: YOU CAN PUT A NOTE HERE
    %
    % **References**
    %
    % [1] Ref1 ...
    %
    % **Example** C=...
    %
    % See also :class:`Map` :class:`Cost`
    
    %%    Copyright (C)
    %     2018 - Biomedical Imaging Group, EPFL
    %	  Masih Nilchian - masih.nilchian@epfl.ch
    %	  Laurène Donati - laurene.donati@epfl.ch
    %     Citation:  Fast Multiresolution Reconstruction for Cryo-EM, Donati et al.
    %
    %     This program is free software: you can redistribute it and/or modify
    %     it under the terms of the GNU General Public License as published by
    %     the Free Software Foundation, either version 3 of the License, or
    %     (at your option) any later version.
    %
    %     This program is distributed in the hope that it will be useful,
    %     but WITHOUT ANY WARRANTY; without even the implied warranty of
    %     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    %     GNU General Public License for more details.
    %
    %     You should have received a copy of the GNU General Public License
    %     along with this program.  If not, see <http://www.gnu.org/licenses/>.
    
    % Protected Set and public Read properties
    properties (SetAccess = protected,GetAccess = public)
        % TODO SET HERE NEW PROTECTED SET AND PUBLIC READ PROPERTIES IF NEEDED.
        
        % PARAMETERS
        param;
        
        % LOOK-UP TABLE
        tau;
        lt;
        gradLT; % the lookup table used for the gradient part
        kernel;
        Ti; % sampling rate of the volume
        Tp; % sampling rate of the projection images
    end
    % Full protected properties
    properties (SetAccess = protected,GetAccess = protected)
        % TODO SET HERE NEW FULLY PROTECTED PROPERTIES
        % (E.G. INTERNAL VARIABLE USED TO AVOID MULTIPLE COMPUTATION)
    end
    
    %% Constructor
    methods
        function this = LinOpPBTShift(sz,angles,shifts,proj_only,kbwf_config) %change the name of the constructor
            
            this.name               ='LinOpPBTShift';
            this.isInvertible       = false;
            this.isDifferentiable   = true;
            this.sizein             = sz;
            
            % KBWF PARAMETERS
            this.param.a               = kbwf_config.a;
            this.param.alpha           = kbwf_config.alpha;
            this.param.m               = kbwf_config.m;
            this.param.scaleRatio      = kbwf_config.scale;
            this.param.aliasingFilter  = false;
            this.param.k1              = this.sizein(1);
            this.param.k2              = this.sizein(2);
            this.param.k3              = this.sizein(3);
            
            % ORGANISE PROJECTION ORIENTATIONS
            this.param.shifts = shifts;
            rotMatrix = eye(3);
            [r1,r2]  = setProjectionPlaneOrientationsSamples(angles, rotMatrix);
            this.param.r1 = r1;
            this.param.r2 = r2;
            
            % CREATE LOOK-UP TABLE
            this.tau = 0.005;
            this.lt = KaiserBesselProjection(this.param.m, this.param.alpha, this.param.a, 0:(this.tau):(this.param.a));

            if proj_only~=1 % if only projections are required and we do not use the operator for reconstruction, then no HTH is required
                this.kernel = fftn(ifftshift(computeHTH_KB_ConventionalCT_3D_r1r2_general(this.param)));
                this.norm = sqrt(max(abs(this.kernel(:))));
            end
            this.Ti = kbwf_config.Ti; this.Tp = kbwf_config.Tp;
            p = this.sizein(1)*this.param.scaleRatio; %sqrt((this.sizein(1)*this.sizein(1)+this.sizein(2)*this.sizein(2)+this.sizein(3)*this.sizein(3))*(scaleRatio^2));
            this.sizeout = [this.sizein(1),this.sizein(1), size(angles,1)];
            if mod(this.sizeout(1), 2)~=0
                this.sizeout(1) = this.sizeout(1)+1;
                this.sizeout(2) = this.sizeout(2)+1;
            end
            if (mod(this.sizeout(1),2)~=0)
                this.sizeout(1) = this.sizeout(1)+1;
                this.sizeout(2) = this.sizeout(2)+1;
            end
        end
    end
    
    %% Core Methods containing implementations (Protected)
    % - apply_(this,x)
    % - applyAdjoint_(this,x)
    methods (Access = protected)
        function y = apply_(this,x)
            % Reimplemented from parent class :class:`Cost`.
            % COMPUTE PROJECTIONS (WITH MULTI-THREAD MEX FILE)
            y = projectionShift3DMT2(x, this.param.r1, this.param.r2, this.param.scaleRatio*[1,1,1]*this.Ti, ...
                [1,1]*this.Tp, this.lt, this.param.a, this.param.shifts(:,1), this.param.shifts(:,2));

        end
        
        function y = applyAdjoint_(this,x)
            % Reimplemented from parent class :class:`Cost`.            
            % COMPUTE HTG
            y = projectionAdjointShift3DMT(x, [this.param.k1; this.param.k2; this.param.k3],...
                this.param.r1, this.param.r2, this.param.scaleRatio*[1,1,1]*this.Ti,...
                [1,1]*this.Tp, this.lt, this.param.a, this.param.shifts(:,1), this.param.shifts(:,2));
        end
        
        function y = applyHtH_(this,x)
            % Reimplemented from parent class :class:`LinOp`.
            y = real(ifftn(fftn(x,[2*this.param.k1-1,2*this.param.k2-1,2*this.param.k3-1]).*this.kernel));
            % crop the region of interest
            y = y(1:this.sizein(1),1:this.sizein(2),1:this.sizein(3));
        end
        
        function M = makeHtH_(this)
            % Reimplemented from parent class :class:`LinOp`.
            S=LinOpSelectorPatch(this.sizein*2-1,[1 1 1],this.sizein);
            H = LinOpConv(this.kernel);
            M = S*H*S';
        end
        
    end
end
