
classdef VM_structure < handle
    %VM_STRUCTURE Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        nex % number of elements in x direction
        ney % number of elements in y direction
        LX % length of mesh domain
        LY % length of mesh domain
        lambda
        mu
        Xmax
        Ymax
        h
        rho0
        numcalls
        error_caught
        num_defects
    end
    
    methods
        function obj = VM_structure(nex,ney,LX,LY,num_defects)
            
            obj.numcalls=0;
            obj.error_caught=0;
            %VM_STRUCTURE Construct an instance of this class
            %   Detailed explanation goes here
            % Parameters for mesh:
            obj.nex = nex; obj.ney = ney; obj.LX = LX; obj.LY = LY;
            obj.num_defects = num_defects;
            Xmin = 0; % minimum X coordinate
            Ymin = 0; % minimum Y coordinate
            obj.Xmax = Xmin+LX; % maximum X coordinate
            obj.Ymax = Ymin+LY; % maximum Y coordinate
            obj.h = LX/nex; % mesh size
            
            % Parameters for material:
            obj.lambda = 100e6; % elastic modulus
            obj.mu = 40e6; % elastic modulus
            obj.rho0 = 0e3; %density
        end

        %parse function for first #defects inputs location of defects.
        function vdata = parse_input(obj, v_input)
            %start parsing first 12
            index_array = 1:25; 
            % prev line: array we are deleting from, the remaining is the 0, deleted is 1
            void_indices = zeros(1,obj.num_defects);
            v25 = zeros(1,25);
            v_def = v_input(1:obj.num_defects); % frist #def elements from our input
            for i = 1: obj.num_defects
                
                void_indices(i) = index_array(v_def(i));
                index_array(v_def(i)) = [];
            end
            v25(void_indices) = 1;
            % end parsing first 12
            vdata = zeros(1,100);
            vdata(1:25) = v25;
            count_defects = 1;
            for i = 1:25
                if v25(i) >= 0.5 %defect in place    
                    vdata(2*i+24 : 2*i+25) = [v_input(obj.num_defects+count_defects) v_input(obj.num_defects+obj.num_defects+count_defects)];
                    vdata(75+i) = v_input(obj.num_defects+2*obj.num_defects+count_defects);
                    count_defects = count_defects + 1;
                else % no defect in place
                    vdata(2*i+24 : 2*i+25) = [0,0];
                    vdata(75+i) = 0;
                end
            end

            %vdata = [v25,v_input(obj.num_defects+1:end)];
        end

        % Computation function
        function cost = compute(obj,vdata)
            try
                displ = -1e-3;
                g = 9.8;
                rho0 = obj.rho0;
                Xmin = 0;
                Ymin = 0;
                
                % v1: ellipse True/False
                % v2: long and short axis length
                % v3: -pi/2 to pi/2
                %vdata = obj.parse(vdata);
                % Display settings
                displayerror = 1; % whether display error message
                displayface = 1; % select face color; 1: interp, 2: texturemap
                displayedge = 0; % select edge color; 0: none, 1: white, 2: black
                hidenodes = 0; % hide nodes near singularities
%% parse vdata
                vdata = obj.parse_input(vdata);
                if sum(size(vdata)) == 2
                    % not satisfying constraint
                    error('vdata not satisfying constraint')
                elseif find_solid_area(vdata) > 0.85
                    disp('solid area exceeding limit')
                end
                v1 = vdata(1:25);
                v2 = vdata(26:75);
                v3 = vdata(76:100);
%% end of parsing
                %vdata = [v1,v2,v3];
                gi = find(v1); % non-zero grid indices
          
                % Elliptical boundaries
%                 ngx = 5; % number of grids in x direction
%                 ngy = 5; % number of grids in y direction
%                 nN = length(gi); % number of Neumann ellipses
%                 XN = Xmin+obj.h*obj.nex/ngx*(floor((gi-1)/ngy)+1/2); % x coordinates of centers of Neumann ellipses
%                 YN = Ymin+obj.h*obj.ney/ngy*(mod(gi-1,ngy)+1/2); % y coordinates of centers of Neumann ellipses
%                 AN = v2(2*gi-1); % semi-major axes of Neumann ellipses
%                 BN = v2(2*gi); % semi-minor axes of Neumann ellipses
%                 thetaN = v3(gi); % inclined angles of Neumann ellipses

                % Basic parameters
%                 nnx = obj.nex+1; % number of nodes in x direction
%                 nny = obj.ney+1; % number of nodes in y direction
%                 ne = obj.nex*obj.ney; % total number of elements
%                 nn = (obj.nex+1)*(obj.ney+1); % total number of nodes
%                 dof = 2*nn; % total degrees of freedom
%                 lcgpe = 1/sqrt(3)*[-1 -1; 1 -1; -1 1; 1 1]'; % Gaussian points in local coordinates
%                 we = 1; % Gaussian weight
                epsilon = 1e-12; % epsilon for crossing points
%                 maxit = 32; % number of maximum iterations
%                 alpha = 1e7; % constant to normalize Dirichlet BC

                if abs(obj.LX/obj.nex-obj.LY/obj.ney) > epsilon && displayerror == 1
                    error('Error in meshing. Elements must be square-shaped.')
                end

                % Elliptical boundaries
                ngx = 5; % number of grids in x direction
                ngy = 5; % number of grids in y direction
                nN = length(gi); % number of Neumann ellipses
                XN = Xmin+obj.h*obj.nex/ngx*(floor((gi-1)/ngy)+1/2); % x coordinates of centers of Neumann ellipses
                YN = Ymin+obj.h*obj.ney/ngy*(mod(gi-1,ngy)+1/2); % y coordinates of centers of Neumann ellipses
                AN = v2(2*gi-1); % semi-major axes of Neumann ellipses
                BN = v2(2*gi); % semi-minor axes of Neumann ellipses
                thetaN = v3(gi); % inclined angles of Neumann ellipses

                % Basic parameters
                nnx = obj.nex+1; % number of nodes in x direction
                nny = obj.ney+1; % number of nodes in y direction
                ne = obj.nex*obj.ney; % total number of elements
                nn = (obj.nex+1)*(obj.ney+1); % total number of nodes
                dof = 2*nn; % total degrees of freedom
                lcgpe = 1/sqrt(3)*[-1 -1; 1 -1; -1 1; 1 1]'; % Gaussian points in local coordinates
                we = 1; % Gaussian weight
                epsilon = 1e-12; % epsilon for crossing points
                maxit = 32; % number of maximum iterations
                alpha = 1e7; % constant to normalize Dirichlet BC

                if abs(obj.LX/obj.nex-obj.LY/obj.ney) > epsilon && displayerror == 1
                    error('Error in meshing. Elements must be square-shaped.')
                end

                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

                % MESHING

                % Connectivity array
                connectivity = zeros(ne,4); % initialize connectivity array
                for n = 1:4 % loop over element nodes
                    for i = 1:obj.nex
                        for j = 1:obj.ney
                            connectivity(i+obj.nex*(j-1),n) = mod(n-1,2)+i+(obj.nex+1)*(floor((n-1)/2)+j-1); % insert node number
                        end
                    end
                end

                % Coordinates array
                coordinates = zeros(2,nn); % initialize coordinates array
                for i = 0:obj.nex
                    for j = 0:obj.ney
                        coordinates(:,1+i+(obj.nex+1)*j) = [Xmin; Ymin] + obj.h*[i;j]; % insert coordinate
                    end
                end

                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

                % LOCAL ARRAYS

                % Interpolation arrays
                Ne = zeros(8,8); % initialize [Ne] at Gaussian points
                Be = zeros(12,8); % initialize [Be] at Gaussian points
                for gp = 1:4 % loop over Gaussian points
                    for n = 1:4 % loop over element nodes
                        xigp = lcgpe(1,gp);
                        etagp = lcgpe(2,gp);
                        Ne(2*gp-1:2*gp,2*n-1:2*n) = 1/4*(1+(-1)^n*xigp)*(1+(-1)^floor((n+1)/2)*etagp)*eye(2);
                        Be(3*gp-2,2*n-1) = 2/obj.h*1/4*(-1)^n*(1+(-1)^floor((n+1)/2)*etagp);
                        Be(3*gp-1,2*n) = 2/obj.h*1/4*(-1)^floor((n+1)/2)*(1+(-1)^n*xigp);
                        Be(3*gp,2*n-1) = 2/obj.h*1/4*(-1)^floor((n+1)/2)*(1+(-1)^n*xigp);
                        Be(3*gp,2*n) = 2/obj.h*1/4*(-1)^n*(1+(-1)^floor((n+1)/2)*etagp);
                    end
                end

                % Element elasticity matrix
                D = [obj.lambda+2*obj.mu obj.lambda 0; obj.lambda obj.lambda+2*obj.mu 0; 0 0 obj.mu]; % [D] (3x3)

                % Element stiffness matrix
                Kefull = zeros(8,8); % initialize [Ke] (8x8)
                for gp = 1:4 % summation over Gaussian points
                    Kefull = Kefull + we*Be(3*gp-2:3*gp,:)'*D*Be(3*gp-2:3*gp,:)*obj.h^2/4;
                end

                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

                % GLOBAL ARRAYS

                % Initialize arrays and other variables
                spi = zeros(64*ne,1); % initialize sparse i
                spj = zeros(64*ne,1); % initialize sparse j
                spv = zeros(64*ne,1); % initialize sparse v
                F = zeros(dof,1); % initialize global forcing vector
                Dfulln = zeros(1,4*ne); % initialize Dirichlet full element numbers
                clon = zeros(1,4*ne); % initialize interior node numbers
                plon = zeros(1,4*ne); % initialize plottable node numbers
                actn = zeros(1,4*ne); % initialize active node numbers
                typen = [0 1 2 3 4; zeros(1,5)]'; % initialize type numbers
                percent  = 0; % initialize progress
                tic

                % Assemble global arrays
                for e = 1:ne % loop over elements

                    % Element information
                    en = connectivity(e,:); % element node numbers
                    edof = sort([2*en-1 2*en]); % element degrees of freedom
                    gcen = coordinates(:,en); % global coordinates of element nodes
                    xmin = gcen(1,1);
                    ymin = gcen(2,1);
                    xmax = gcen(1,4);
                    ymax = gcen(2,4);
                    xc = (xmin+xmax)/2;
                    yc = (ymin+ymax)/2;
                    gcgp = gcen*Ne([1 3 5 7],[1 3 5 7])'; % global coordinates of Gaussian points

                    % Interior and exterior nodes
                    nint = 0; % initialize number of interior nodes
                    next = 0; % initialize number of exterior nodes
                    intn = zeros(1,4); % initialize interior node numbers
                    extn = zeros(1,4); % initialize exterior node numbers
                    Nn = 0; % initialize Neumann boundary number
                    for n = 1:4 % loop over element nodes
                        xn = gcen(1,n); % x coordinate of node
                        yn = gcen(2,n); % y coordinate of node
                        if xn^2+yn^2 >= 0 % input solid boundary inequality (interior)
                            Nntemp = Nn;
                            for i = 1:nN % loop over Neumann boundaries
                                if ((xn-XN(i))*cos(thetaN(i))+(yn-YN(i))*sin(thetaN(i)))^2/AN(i)^2+(-(xn-XN(i))*sin(thetaN(i))+(yn-YN(i))*cos(thetaN(i)))^2/BN(i)^2-1 < 0 % input Neumann hollow boundary inequality (exterior)
                                    if Nntemp ~= 0 && Nntemp ~= nN+1 && Nntemp ~= i && displayerror == 1 % element invloves more than one Neumann hollow boundaries
                                        error(['Error in classifying partial node ',num2str(en(n)),' at ',mat2str(gcen(:,n)'),'.'])
                                    end
                                    Nn = i; % Neumann boundary number
                                    Nntemp = Nn + 1;
                                end
                            end
                            if Nn == Nntemp % interior node
                                nint = nint + 1;
                                intn(nint) = en(n); % interior node number
                            elseif Nn == Nntemp-1 % Neumann exterior node
                                next = next + 1;
                                extn(next) = en(n); % exterior node number
                            elseif displayerror == 1
                                error(['Error in classifying partial node ',num2str(en(n)),' at ',mat2str(gcen(:,n)'),'.'])
                            end
                        else
                            next = next + 1;
                            extn(next) = en(n); % exterior node number
                        end
                    end
                    intn = sort(intn(1:nint)); % interior node numbers (ascending order)
                    extn = sort(extn(1:next)); % exterior node numbers (ascending order)

                    % Element arrays
                    Fbefull = zeros(8,1); % initialize full element body force vector
                    Ftefull = zeros(8,1); % initialize full element traction vector
                    if nint > 0 % active node
                        for gp = 1:4 % summation over Gaussian points
%                             xgp = gcgp(1,gp); % x coordinate of Gaussian point
%                             ygp = gcgp(2,gp); % y coordinate of Gaussian point
                            bgp = [0; -obj.rho0*g]; % input body force
                            if norm(bgp) > 0
                                Fbefull = Fbefull + we*Ne(2*gp-1:2*gp,:)'*bgp*h^2/4;
                            end
                        end
                        if xmin == Xmin % input full element Neumann boundary equation
                            xgp = xmin; % x coordinate of Gaussian point
                            Nb = zeros(4,8); % initialize [Nb] at Gaussian points
                            for gp = 1:2 % loop over Gaussian points
                                ygp = gcgp(2,2*gp-1); % y coordinate of Gaussian point
                                for n = 1:4 % loop over element nodes
                                    xigp = (xgp-xc)*2/obj.h;
                                    etagp = (ygp-yc)*2/obj.h;
                                    Nb(2*gp-1:2*gp,2*n-1:2*n) = 1/4*(1+(-1)^n*xigp)*(1+(-1)^floor((n+1)/2)*etagp)*eye(2);
                                end
                                tgp = [0; 0]; % input traction
                                if norm(tgp) > 0
                                    Ftefull = Ftefull + we*Nb(2*gp-1:2*gp,:)'*tgp*obj.h/2;
                                end
                            end
                        end
                        if xmax == obj.Xmax % input full element Neumann boundary equation
                            xgp = xmax; % x coordinate of Gaussian point
                            Nb = zeros(4,8); % initialize [Nb] at Gaussian points
                            for gp = 1:2 % loop over Gaussian points
                                ygp = gcgp(2,2*gp-1); % y coordinate of Gaussian point
                                for n = 1:4 % loop over element nodes
                                    xigp = (xgp-xc)*2/obj.h;
                                    etagp = (ygp-yc)*2/obj.h;
                                    Nb(2*gp-1:2*gp,2*n-1:2*n) = 1/4*(1+(-1)^n*xigp)*(1+(-1)^floor((n+1)/2)*etagp)*eye(2);
                                end
                                tgp = [0; 0]; % input traction
                                if norm(tgp) > 0
                                    Ftefull = Ftefull + we*Nb(2*gp-1:2*gp,:)'*tgp*obj.h/2;
                                end
                            end
                        end
                        if ymin == Ymin % input full element Neumann boundary equation
                            ygp = ymin;
                            Nb = zeros(4,8); % initialize [Nb] at Gaussian points
                            for gp = 1:2 % loop over Gaussian points
                                xgp = gcgp(1,gp); % x coordinate of Gaussian point
                                for n = 1:4 % loop over element nodes
                                    xigp = (xgp-xc)*2/obj.h;
                                    etagp = (ygp-yc)*2/obj.h;
                                    Nb(2*gp-1:2*gp,2*n-1:2*n) = 1/4*(1+(-1)^n*xigp)*(1+(-1)^floor((n+1)/2)*etagp)*eye(2);
                                end
                                tgp = [0; 0]; % input traction
                                if norm(tgp) > 0
                                    Ftefull = Ftefull + we*Nb(2*gp-1:2*gp,:)'*tgp*obj.h/2;
                                end
                            end
                        end
                        if ymax == obj.Ymax % input full element Neumann boundary equation
                            ygp = ymax;
                            Nb = zeros(4,8); % initialize [Nb] at Gaussian points
                            for gp = 1:2 % loop over Gaussian points
                                xgp = gcgp(1,gp); % x coordinate of Gaussian point
                                for n = 1:4 % loop over element nodes
                                    xigp = (xgp-xc)*2/obj.h;
                                    etagp = (ygp-yc)*2/obj.h;
                                    Nb(2*gp-1:2*gp,2*n-1:2*n) = 1/4*(1+(-1)^n*xigp)*(1+(-1)^floor((n+1)/2)*etagp)*eye(2);
                                end
                                tgp = [0; 0]; % input traction
                                if norm(tgp) > 0
                                    Ftefull = Ftefull + we*Nb(2*gp-1:2*gp,:)'*tgp*obj.h/2;
                                end
                            end
                        end
                    end
                    clon(4*e-3:4*(e-1)+nint) = intn; % interior node numbers
                    if nint == 0 % empty element
                        typen(1,2) = typen(1,2) + 1;
                        Ke = zeros(8,8); % empty element stiffness matrix
                        Fbe = zeros(8,1); % empty element body force vector
                        Fte = zeros(8,1); % empty element traction vector
                    elseif nint == 4 % full element
                        typen(5,2) = typen(5,2) + 1;
                        Ke = Kefull; % full element stiffness matrix
                        Fbe = Fbefull; % full element body force vector
                        Fte = Ftefull; % full element traction vector
                        for n = 1:4 % loop over element nodes
                            if abs(gcen(1,n)-Xmin) < epsilon || abs(gcen(1,n)-obj.Xmax) < epsilon || abs(gcen(2,n)-Ymin) < epsilon || abs(gcen(2,n)-obj.Ymax) < epsilon % input full element Dirichlet boundary equation
                            %if gcen(1,n) == Xmin || gcen(1,n) == Xmax || gcen(2,n) == Ymin || gcen(2,n) == Ymax % input full element Dirichlet boundary equation
                                Dfulln(4*(e-1)+n) = en(n);
                            end
                        end
                        plon(4*e-3:4*e) = en; % plottable node numbers
                        actn(4*e-3:4*e) = en; % active node numbers
                    else % partial element
                        gcint = coordinates(:,intn); % global coordinates of interior nodes (ascending order)
                        gcext = coordinates(:,extn); % global coordinates of exterior nodes (ascending order)
                        if Nn > 0 % Neumann partial element
                            [Ke,Fbe,Fte,typen] = NeumannBC(obj,XN(Nn),YN(Nn),AN(Nn),BN(Nn),thetaN(Nn),obj.h,e,en,gcen,xmin,ymin,xmax,ymax,nint,next,intn,extn,typen,D,obj.rho0,g,Kefull,Fbefull,maxit,epsilon,displayerror); % Neumann partial element arrays
                            plon(4*e-3:4*(e-1)+nint) = intn; % plottable node numbers
                        elseif displayerror == 1
                            error(['Error in classifying partial element ',num2str(e),' between ',mat2str(gcen(:,1)'),' and ',mat2str(gcen(:,4)'),'.'])
                        end
                        actn(4*e-3:4*e) = en; % active node numbers
                    end

                    % Assemble global arrays
                    spi(64*e-63:64*e) = repmat(edof',[8 1]); % insert element row indices
                    spj(64*e-63:64*e) = repelem(edof',8); % insert element column indices
                    spv(64*e-63:64*e) = Ke(:); % insert element stiffness values
                    F(edof) = F(edof) + Fbe + Fte; % insert element forcing vectors

                    % Show progress
                    if floor(100*e/ne) > percent
                        percent = floor(100*e/ne);
                    end

                end

                K = sparse(spi,spj,spv,dof,dof); % form global stiffness matrix
                K0 = K; % global stiffness matrix prior to imposing Dirichlet BC


                % Impose Dirichlet BC on full elements
                Dfulln = setdiff(Dfulln,0); % Dirichlet full element nodes
                Dfulldof = zeros(2*length(Dfulln),1);
                Dfullubar = zeros(2*length(Dfulln),1);
                if ~isempty(Dfulln) % Dirichlet full element nodes are non-empty
                    for i = 1:length(Dfulln) % loop over Dirichlet full element nodes
                        n = Dfulln(i);
                        xn = coordinates(1,Dfulln(i));
                        yn = coordinates(2,Dfulln(i));
                        if n == nnx
                            Dc = 1:2; % input Dirichlet BC coordinates
                            ubar = [0; 0]; % input Dirichlet BC
                        elseif n == nnx*nny
                            Dc = 1:2; % input Dirichlet BC coordinates
                            ubar = [0; displ]; % input Dirichlet BC
                        elseif n <= nnx
                            Dc = 2; % input Dirichlet BC coordinates
                            ubar = [0; 0]; % input Dirichlet BC
                        elseif n >= nnx*(nny-1)+1
                            Dc = 2; % input Dirichlet BC coordinates
                            ubar = [0; displ]; % input Dirichlet BC
                        elseif mod(n,nnx) == 0
                            Dc = 1; % input Dirichlet BC coordinates
                            ubar = [0; 0]; % input Dirichlet BC
                        else
                            Dc = 0; % no Dirichlet BC
                        end
                        if find(Dc==1)
                            Dfulldof(2*i-1) = 2*Dfulln(i)-1;
                            Dfullubar(2*i-1) = ubar(1);
                        end
                        if find(Dc==2)
                            Dfulldof(2*i) = 2*Dfulln(i);
                            Dfullubar(2*i) = ubar(2);
                        end
                    end
                    Dfulldofi = find(Dfulldof~=0);
                    Dfulldof = Dfulldof(Dfulldofi);
                    Dfullubar = Dfullubar(Dfulldofi);
                    K(Dfulldof,:) = sparse(length(Dfulldof),dof); % initialize rows of global stiffness matrix associated with Dirichlet full element nodes
                    K(Dfulldof,Dfulldof) = alpha*speye(length(Dfulldof)); % insert identity matrix
                    F(Dfulldof) = alpha*Dfullubar; % insert Dirichlet BC
                end

                % Types of nodes
                clon = setdiff(clon,0); % interior node numbers
                comn = setdiff(1:nn,clon); % exterior node numbers
                plon = setdiff(plon,0); % plottable node numbers
                unpn = setdiff(1:nn,plon); % hidden node numbers
                actn = setdiff(actn,0); % active node numbers
                inan = setdiff(1:nn,actn); % inactive node numbers
                actdof = sort([2*actn-1 2*actn]); % active degress of freedom
                Kact = K(actdof,actdof); % K containing active dofs only
                K0act = K0(actdof,actdof);
                Fact = F(actdof); % F containing active dofs only

                Kact = sparse(Kact);
                K0act = sparse(K0act);
                Fact = sparse(Fact);

                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

                % SOLUTIONS AND PLOTS


                % Solutions
                uact = Kact\Fact; % solve for active nodes

                toc
                % Process solutions
                u = zeros(dof,1);
                u(actdof) = uact;
                u(sort([2*comn-1 2*comn])) = NaN;
                r = zeros(dof,1);
                r(actdof) = K0act*uact;
                r(sort([2*comn-1 2*comn])) = NaN;

                % Displacement and stress data
                u1u2 = reshape(u,[2,nn]);
                u1 = u1u2(1,:)';
                u2 = u1u2(2,:)';
                sn = sqrt(u1.^2+u2.^2);
                r1r2 = reshape(r,[2,nn]);
                r1 = r1r2(1,:)';
                r2 = r1r2(2,:)';
                epsilon11 = zeros(ne,1);
                epsilon22 = zeros(ne,1);
                epsilon33 = zeros(ne,1);
                epsilon12 = zeros(ne,1);
                epsilon23 = zeros(ne,1);
                epsilon31 = zeros(ne,1);
                sigma11 = zeros(ne,1);
                sigma22 = zeros(ne,1);
                sigma33 = zeros(ne,1);
                sigma12 = zeros(ne,1);
                sigma23 = zeros(ne,1);
                sigma31 = zeros(ne,1);
                sigmavM = zeros(ne,1);
                for e = 1:ne % loop over elements
                    en = connectivity(e,:); % element node numbers
                    epsilon11(e) = (u1(en(2))+u1(en(4))-u1(en(1))-u1(en(3)))/(2*obj.h);
                    epsilon22(e) = (u2(en(3))+u2(en(4))-u2(en(1))-u2(en(2)))/(2*obj.h);
                    epsilon12(e) = (u1(en(3))+u1(en(4))-u1(en(1))-u1(en(2))+u2(en(2))+u2(en(4))-u2(en(1))-u2(en(3)))/(4*obj.h);
                    sigma11(e) = obj.lambda*(epsilon11(e)+epsilon22(e)+epsilon33(e))+2*obj.mu*epsilon11(e);
                    sigma22(e) = obj.lambda*(epsilon11(e)+epsilon22(e)+epsilon33(e))+2*obj.mu*epsilon22(e);
                    sigma33(e) = obj.lambda*(epsilon11(e)+epsilon22(e)+epsilon33(e))+2*obj.mu*epsilon33(e);
                    sigma12(e) = 2*obj.mu*epsilon12(e);
                    sigma23(e) = 2*obj.mu*epsilon23(e);
                    sigma31(e) = 2*obj.mu*epsilon31(e);
                    sigmavM(e) = sqrt(((sigma11(e)-sigma22(e))^2+(sigma22(e)-sigma33(e))^2+(sigma33(e)-sigma11(e))^2+6*(sigma12(e)^2+sigma23(e)^2+sigma31(e)^2))/2);
                end
                % Reaction nodal forces
                F1bottom = sum(r1(1:nnx));
                F1total = sum(r1,'omitnan');
                F2bottom = sum(r2(1:nnx));
                F2top = sum(r2(nnx*(nny-1)+1:nnx*nny));
                F2row = sum(r2(nnx+1:2*nnx));
                F2total = sum(r2,'omitnan');
                % Reaction surface forces
                disp('Reaction forces:')
                R2top = F2top-F2row/2;
                R2bottom = F2bottom-F2row/2;
                FEMweight=R2top+R2bottom;
                % Weight
                Vrect = obj.LX*obj.LY;
                Vvoid = pi*dot(AN,BN);
                Vmatl = Vrect-Vvoid;
                weight = obj.rho0*g*Vmatl;
                % Stiffness
                stiffnesspervolume = (abs(R2top)+abs(R2bottom))/2/abs(displ)/Vmatl;
                % von Mises stress
                nactn = obj.nex*obj.ney-length(find(isnan(sigmavM)));
                avgvM = sum(sigmavM,'omitnan')/nactn;
                maxvM = max(sigmavM);
                minvM = min(sigmavM);
                %disp('Finished solving.')
% 
%                 % Process solutions
%                 u = zeros(dof,1);
%                 u(actdof) = uact;
%                 u(sort([2*comn-1 2*comn])) = NaN;
%                 % Displacement and stress data
%                 u1u2 = reshape(u,[2,nn]);
%                 u1 = u1u2(1,:)';
%                 u2 = u1u2(2,:)';
%                 sn = sqrt(u1.^2+u2.^2);
%                 epsilon11 = zeros(ne,1);
%                 epsilon22 = zeros(ne,1);
%                 epsilon33 = zeros(ne,1);
%                 epsilon12 = zeros(ne,1);
%                 epsilon23 = zeros(ne,1);
%                 epsilon31 = zeros(ne,1);
%                 sigma11 = zeros(ne,1);
%                 sigma22 = zeros(ne,1);
%                 sigma33 = zeros(ne,1);
%                 sigma12 = zeros(ne,1);
%                 sigma23 = zeros(ne,1);
%                 sigma31 = zeros(ne,1);
%                 sigmavM = zeros(ne,1);
%                 for e = 1:ne % loop over elements
%                     en = connectivity(e,:); % element node numbers
%                     epsilon11(e) = (u1(en(2))+u1(en(4))-u1(en(1))-u1(en(3)))/(2*obj.h);
%                     epsilon22(e) = (u2(en(3))+u2(en(4))-u2(en(1))-u2(en(2)))/(2*obj.h);
%                     epsilon12(e) = (u1(en(3))+u1(en(4))-u1(en(1))-u1(en(2))+u2(en(2))+u2(en(4))-u2(en(1))-u2(en(3)))/(4*obj.h);
%                     sigma11(e) = obj.lambda*(epsilon11(e)+epsilon22(e)+epsilon33(e))+2*obj.mu*epsilon11(e);
%                     sigma22(e) = obj.lambda*(epsilon11(e)+epsilon22(e)+epsilon33(e))+2*obj.mu*epsilon22(e);
%                     sigma33(e) = obj.lambda*(epsilon11(e)+epsilon22(e)+epsilon33(e))+2*obj.mu*epsilon33(e);
%                     sigma12(e) = 2*obj.mu*epsilon12(e);
%                     sigma23(e) = 2*obj.mu*epsilon23(e);
%                     sigma31(e) = 2*obj.mu*epsilon31(e);
%                     sigmavM(e) = sqrt(((sigma11(e)-sigma22(e))^2+(sigma22(e)-sigma33(e))^2+(sigma33(e)-sigma11(e))^2+6*(sigma12(e)^2+sigma23(e)^2+sigma31(e)^2))/2);
%                 end
%                 nactn = obj.nex*obj.ney-length(find(isnan(sigmavM)));
%                 avgvM = sum(sigmavM,'omitnan')/nactn;
%                 concrete_area_percentage = nactn/(obj.nex*obj.ney);
%                 %disp([num2str(concrete_area_percentage*100),'% concrete area'])
%                 
%                 %edit 8/24/22 for centroid
%                 sigmavMc = reshape(sigmavM,obj.nex,obj.ney)';
%                 geometry_matrix = (isnan(sigmavMc)*(-1))+1;
%                 centroids = regionprops(true(size(geometry_matrix)), geometry_matrix, 'WeightedCentroid');
%                 x_mean_node = centroids.WeightedCentroid(1);
%                 y_mean_node = centroids.WeightedCentroid(2);
%                 %disp(['y_mean_node: ', num2str(y_mean_node)])
%                 %end of edit 8/24/22
%                 
%                 %edit 8/25/22 for strain displacement
%           
%                 index_of_centroid = sub2ind([obj.nex, obj.ney], obj.nex, floor(y_mean_node));
%                 centroid_strain = epsilon22(index_of_centroid);
%                 %disp(['strain at center is :',centroid_strain])
%                 %end of edit 8/25/22
%                 
%                 % edit for cost definition and opt
%                 temp = sum(sigmavM,'omitnan')/3e10;
%                 %cost = norm([concrete_area_percentage, temp]);
%                 cost = abs(centroid_strain); %cost edit 8/25/22
%                 %disp(['sigma % is ', num2str(temp*100)])
%                 %maxvM = max(sigmavM);
%                 %minvM = min(sigmavM);
%                 %disp(['Average von Mises stress = ',num2str(avgvM),' Pa'])
%                 %disp(['Cost of Opt = ',num2str(cost)])
%                 VonMisesStress = avgvM;



                %for recording only
                cost = abs(stiffnesspervolume)/1e7; %normalized to up to 3.3
                disp(['cost is : ',num2str(cost)])
                vdata = [v1,v2,v3];
                writematrix([vdata'; cost],['sim_results/sim',num2str(obj.numcalls)])
                disp(['results in sim',num2str(obj.numcalls)])
                obj.change_numcalls(obj.numcalls + 1);
                %end of edit for cost definition and opt
            catch ME
                
                obj.change_error_caught(obj.error_caught + 1);
                errorMessage = sprintf('Error in function %s() at line %d.\n\nError Message:\n%s', ...
                ME.stack(1).name, ME.stack(1).line, ME.message);
                fprintf(1, '%s\n', errorMessage);
                cost = 0; %for unexpected errors
                %writematrix([vdata'],['error',num2str(obj.error_caught)])
            end
        end    
        
        
        function vdata = parse(obj, vdata75)
            cur = obj;
            vdata = [zeros(1,25),vdata75];
            for i = 1:2:50
                if vdata75(i)~=0 && vdata75(i+1) ~= 0
                    vdata((i+1)/2) = 1;
                else
                    vdata((i+1)/2) = 0;
                end
            end       
        end
        
        function numcalls = change_numcalls(obj,val)
            obj.numcalls = val;
            numcalls = obj.numcalls;
        end
        
        function error_caught = change_error_caught(obj,val)
            obj.error_caught = val;
            error_caught = obj.error_caught; 
        end
        
        % Boundary Condition Function used in compute
        function [Ke,Fbe,Fte,typen] = NeumannBC(obj,XN,YN,AN,BN,thetaN,h,e,en,gcen,xmin,ymin,xmax,ymax,nint,next,intn,extn,typen,D,rho0,g,Kefull,Fbefull,maxit,epsilon,displayerror)

            % Obtain Neumann crossing points
            gcmin = gcen(:,1);
            gcmax = gcen(:,4);
            gcc = 1/2*(gcmin+gcmax);
            xc = gcc(1);
            yc = gcc(2);
            ncro = 0;
            gccro = zeros(2,4);

            % Newton-Raphson method (ellipse)
            % Solve for x
            sol = [ymin ymax];
            for i = 1:2 % loop over edges
                ysol = sol(i);
                xsol = xc;
                for k = 1:maxit % N-R iterations
                    fx = ((xsol-XN)*cos(thetaN)+(ysol-YN)*sin(thetaN))^2/AN^2+(-(xsol-XN)*sin(thetaN)+(ysol-YN)*cos(thetaN))^2/BN^2-1; % ellipse
                    residue = abs(fx);
                    if residue < epsilon
                        if xmin < xsol && xsol < xmax
                            ncro = ncro + 1;
                            gccro(:,ncro) = [xsol; ysol];
                        end
                        %k
                        break
                    end
                    Dfx = 2*((xsol-XN)*cos(thetaN)+(ysol-YN)*sin(thetaN))*cos(thetaN)/AN^2-2*(-(xsol-XN)*sin(thetaN)+(ysol-YN)*cos(thetaN))*sin(thetaN)/BN^2; % ellipse
                    dx = -fx/Dfx;
                    xsol = xsol + dx;
                    if k == maxit
                        break
                    end
                end
            end
            ncrox = ncro;
            % Solve for y
            sol = [xmin xmax];
            for i = 1:2 % loop over edges
                xsol = sol(i);
                ysol = yc;
                for k = 1:maxit % N-R iterations
                    fy = ((xsol-XN)*cos(thetaN)+(ysol-YN)*sin(thetaN))^2/AN^2+(-(xsol-XN)*sin(thetaN)+(ysol-YN)*cos(thetaN))^2/BN^2-1; % ellipse
                    residue = abs(fy);
                    if residue < epsilon
                        if ymin < ysol && ysol < ymax
                            ncro = ncro + 1;
                            gccro(:,ncro) = [xsol; ysol];
                        end
                        %k
                        break
                    end
                    Dfy = 2*((xsol-XN)*cos(thetaN)+(ysol-YN)*sin(thetaN))*sin(thetaN)/AN^2+2*(-(xsol-XN)*sin(thetaN)+(ysol-YN)*cos(thetaN))*cos(thetaN)/BN^2; % ellipse
                    dy = -fy/Dfy;
                    ysol = ysol + dy;
                    if k == maxit
                        break
                    end
                end
            end
            ncroy = ncro-ncrox;

            if ncro ~= 2
                if nint == 2 && ncrox == 2 && ncroy == 1
                    ncro = 2;
                    gccro = gccro(:,1:2);
                elseif nint == 2 && ncrox == 1 && ncroy == 2
                    ncro = 2;
                    gccro = gccro(:,2:3);
                elseif displayerror == 1
                    error(['Error in identifying crossing points of partial element ',num2str(e),' between ',mat2str(gcen(:,1)'),' and ',mat2str(gcen(:,4)'),'.'])
                end
            end

            gccro = gccro(:,1:ncro); % global coordinates of crossing points (not ordered)

            % Obtain interior or exterior node indices and global coordinates (ascending order)
            if nint <= 2
                inti = zeros(1,nint); % initialize interior node indices
                for i = 1:nint % loop over interior nodes
                    inti(i) = find(en==intn(i)); % insert interior node index
                end
                gcint = gcen(:,inti); % global coordinates of interior nodes
            else
                exti = zeros(1,next); % initialize exterior node indices
                for i = 1:next % loop over exterior nodes
                    exti(i) = find(en==extn(i)); % insert exterior node index
                end
                gcext = gcen(:,exti); % global coordinates of exterior nodes
            end
            gccrotemp = zeros(4,ncro); % initialize temporary matrix

            % Identify cases based on number of interior nodes
            switch nint
                case 1 % one interior node
                    nver = 3; % number of vertices
                    ntri = 1; % number of triangles
                    nseg = 1; % number of segments
                    for i = 1:ncro % loop over Neumann crossing points
                        for j = 1:2 % loop to associate with interior nodes
                            if gccro(mod(j,2)+1,i) == gcint(mod(j,2)+1,1) % Neumann crossing point and interior node are on the same j axis
                                gccrotemp(1:2,j) = gccro(:,i); % insert global coordinates of Neumann crossing points
                                gccrotemp(3,j) = gccrotemp(3,j)+1;
                            end
                        end
                    end
                    if ~isequal(gccrotemp(3,:),ones(size(gccrotemp(3,:)))) && displayerror == 1
                        error(['Error in identifying Neumann crossing points of partial element ',num2str(e),' between ',mat2str(gcen(:,1)'),' and ',mat2str(gcen(:,4)'),'.'])
                    end
                    gccro = gccrotemp(1:2,:); % global coordinates of Neumann crossing points (re-ordered)
                    gcpol = zeros(2,nver); % initialize global coordinates of triangular vertices
                    gcpol(:,1) = gcint; % insert global coordinates of interior node
                    gcpol(:,2:3) = gccro; % insert global coordinates of Neumann boundary
                    trii = [1 2 3]; % triangular vertex indices
                    segi = [2 3]; % segmental vertex indices
                    typen(2,2) = typen(2,2) + 1;
                case 2 % two interior nodes
                    nver = 4; % number of vertices
                    ntri = 2; % number of triangles
                    nseg = 1; % number of segments
                    for i = 1:ncro % loop over Neumann crossing points
                        for j = 1:2 % loop to associate with interior nodes
                            if gccro(mod(j,2)+1,i) == gcint(mod(j,2)+1,1) % Neumann crossing point and interior node are on the same j axis
                                gccrotemp(1:2,1) = gccro(:,i); % insert global coordinates of Neumann crossing points
                                gccrotemp(3,1) = gccrotemp(3,1)+1;
                            elseif gccro(mod(j,2)+1,i) == gcint(mod(j,2)+1,2) % Neumann crossing point and interior node are on the same j axis
                                gccrotemp(1:2,2) = gccro(:,i); % insert global coordinates of Neumann crossing points
                                gccrotemp(3,2) = gccrotemp(3,2)+1;
                            end
                        end
                    end
                    if ~isequal(gccrotemp(3,:),ones(size(gccrotemp(3,:)))) && displayerror == 1
                        error(['Error in identifying Neumann crossing points of partial element ',num2str(e),' between ',mat2str(gcen(:,1)'),' and ',mat2str(gcen(:,4)'),'.'])
                    end
                    gccro = gccrotemp(1:2,:); % global coordinates of Neumann crossing points (re-ordered)
                    gcpol = zeros(2,nver); % initialize global coordinates of triangular vertices
                    gcpol(:,1:2) = gcint; % insert global coordinates of interior node
                    gcpol(:,3:4) = gccro; % insert global coordinates of Neumann boundary
                    trii = [[1 2 3] [2 3 4]]; % triangular vertex indices
                    segi = [3 4]; % segmental vertex indices
                    typen(3,2) = typen(3,2) + 1;
                case 3 % one exterior node
                    nver = 3; % number of vertices
                    ntri = 1; % number of triangles
                    nseg = 1; % number of segments
                    for i = 1:ncro % loop over Neumann crossing points
                        for j = 1:2 % loop to associate with interior nodes
                            if gccro(mod(j,2)+1,i) == gcext(mod(j,2)+1,1) % Neumann crossing point and exterior node are on the same j axis
                                gccrotemp(1:2,j) = gccro(:,i); % insert global coordinates of Neumann crossing points
                                gccrotemp(3,j) = gccrotemp(3,j)+1;
                            end
                        end
                    end
                    if ~isequal(gccrotemp(3,:),ones(size(gccrotemp(3,:)))) && displayerror == 1
                        error(['Error in identifying Neumann crossing points of partial element ',num2str(e),' between ',mat2str(gcen(:,1)'),' and ',mat2str(gcen(:,4)'),'.'])
                    end
                    gccro = gccrotemp(1:2,:); % global coordinates of Neumann crossing points (re-ordered)
                    gcpol = zeros(2,nver); % initialize global coordinates of triangular vertices
                    gcpol(:,1) = gcext; % insert global coordinates of exterior node
                    gcpol(:,2:3) = gccro; % insert global coordinates of Neumann boundary
                    trii = [1 2 3]; % triangular vertex indices
                    segi = [2 3]; % segmental vertex indices
                    typen(4,2) = typen(4,2) + 1;
                otherwise
                    if displayerror == 1
                        error(['Error in classifying partial element ',num2str(e),' between ',mat2str(gcen(:,1)'),' and ',mat2str(gcen(:,4)'),'.'])
                    end
            end
            lcpol = (gcpol-gcc)*2/h; % local coordinates of polygonal vertices

            % Obtain Gaussian points for polygon
            gcpolgp = zeros(2,3*ntri); % initialize global coordinates of Gaussian points of polygon
            for i = 1:ntri % loop over triangles
                for gp = 1:3 % loop over Gaussian points of triangle
                    gcpolgp(:,3*(i-1)+gp) = sum(gcpol(:,trii(setdiff(3*i-2:3*i,3*(i-1)+gp))),2)/2; % insert global coordinates of Gaussian point
                end
            end
            lcpolgp = (gcpolgp-gcc)*2/h; % local coordinates of Gaussian points
            Ntri = zeros(6*ntri,8); % initialize [Ne] at Gaussian points
            Btri = zeros(9*ntri,8); % initialize [Be] at Gaussian points
            for gp = 1:3*ntri % loop over Gaussian points
                xigp = lcpolgp(1,gp);
                etagp = lcpolgp(2,gp);
                for n = 1:4 % loop over element nodes
                    Ntri(2*gp-1:2*gp,2*n-1:2*n) = 1/4*(1+(-1)^n*xigp)*(1+(-1)^floor((n+1)/2)*etagp)*eye(2);
                    Btri(3*gp-2,2*n-1) = 2/h*1/4*(-1)^n*(1+(-1)^floor((n+1)/2)*etagp);
                    Btri(3*gp-1,2*n) = 2/h*1/4*(-1)^floor((n+1)/2)*(1+(-1)^n*xigp);
                    Btri(3*gp,2*n-1) = 2/h*1/4*(-1)^floor((n+1)/2)*(1+(-1)^n*xigp);
                    Btri(3*gp,2*n) = 2/h*1/4*(-1)^n*(1+(-1)^floor((n+1)/2)*etagp);
                end
            end
            wtri = 1/3;
            lctri = lcpol(:,trii); % local coordinates of polygonal vertices
            latri = zeros(1,ntri); % initialize areas of triangles using local coordinates
            for i = 1:ntri % loop over triangles
                latri(i) = 1/2*norm(cross([lctri(:,3*(i-1)+2)-lctri(:,3*(i-1)+1); 0],[lctri(:,3*(i-1)+3)-lctri(:,3*(i-1)+1); 0])); % insert area of triangle using local coordinates
            end

            % Obtain Gaussian points for segment
            gcseggp = zeros(2,2*nseg); % initialize global coordinates of Gaussian points of segment
            for i = 1:nseg % loop over segments
                for gp = 1:2 % loop over Gaussian points of segment
                    gcseggp(:,2*(i-1)+gp) = 1/2*(1-(-1)^gp/sqrt(3))*gcpol(:,segi(2*i-1))+1/2*(1+(-1)^gp/sqrt(3))*gcpol(:,segi(2*i)); % insert global coordinates of Gaussian point
                end
            end
            lcseggp = (gcseggp-gcc)*2/h; % local coordinates of Gaussian points
            Nseg = zeros(4*nseg,8); % initialize [Ne] at Gaussian points
            Bseg = zeros(6*nseg,8); % initialize [Be] at Gaussian points
            for gp = 1:2*nseg % loop over Gaussian points
                xigp = lcseggp(1,gp);
                etagp = lcseggp(2,gp);
                for n = 1:4 % loop over element nodes
                    Nseg(2*gp-1:2*gp,2*n-1:2*n) = 1/4*(1+(-1)^n*xigp)*(1+(-1)^floor((n+1)/2)*etagp)*eye(2);
                    Bseg(3*gp-2,2*n-1) = 2/h*1/4*(-1)^n*(1+(-1)^floor((n+1)/2)*etagp);
                    Bseg(3*gp-1,2*n) = 2/h*1/4*(-1)^floor((n+1)/2)*(1+(-1)^n*xigp);
                    Bseg(3*gp,2*n-1) = 2/h*1/4*(-1)^floor((n+1)/2)*(1+(-1)^n*xigp);
                    Bseg(3*gp,2*n) = 2/h*1/4*(-1)^n*(1+(-1)^floor((n+1)/2)*etagp);
                end
            end
            wseg = 1/2;
            lcseg = lcpol(:,segi); % local coordinates of segmental vertices
            llseg = zeros(1,nseg); % initialize lengths of segments using local coordinates
            for i = 1:nseg % loop over segments
                llseg(i) = norm(lcseg(:,2*i)-lcseg(:,2*i-1)); % insert length of segment using local coordinates
            end

            % Normal traction / pressure
            gcsegc = zeros(2,nseg); % initialize global coordinates of segmental centers
            gcpi = zeros(2,nseg); % initialize global coordinates of pi vectors
            tang = zeros(2,nseg); % initialize tangent vectors
            onor = zeros(2,nseg); % initialize outward normal vectors
            nort = zeros(2,nseg); % initialize normal tractions
            for i = 1:nseg % loop over segments
                gcsegc(:,i) = 1/2*(gcpol(:,segi(2*i-1))+gcpol(:,segi(2*i)));
                gcpi(:,i) = gcsegc(:,i)-[XN; YN];
                tang(:,i) = (gcpol(:,segi(2*i))-gcpol(:,segi(2*i-1)))/norm(gcpol(:,segi(2*i))-gcpol(:,segi(2*i-1)));
                onor(:,i) = [tang(2,i); -tang(1,i)];
                if dot(gcpi,onor) < 0
                    onor(:,i) = -onor(:,i);
                end
                nort(:,i) = onor(:,i)*llseg(i)/2;
            end

            % Obtain local arrays
            if nint <= 2 % evaluate using interior polygon
                Keint = zeros(8,8); % initialize element stiffness matrix for interior domain of partial element
                Fbeint = zeros(8,1); % initialize element body force vector for interior domain of partial element
                for i = 1:ntri % summation over triangles
                    for gp = 1:3 % summation over Gaussian points of triangles
                        Keint = Keint + wtri*Btri(9*(i-1)+3*gp-2:9*(i-1)+3*gp,:)'*D*Btri(9*(i-1)+3*gp-2:9*(i-1)+3*gp,:)*latri(i)*h^2/4;
                        %xgp = gcpolgp(1,3*(i-1)+gp); % x coordinate of Gaussian point
                        %ygp = gcpolgp(2,3*(i-1)+gp); % y coordinate of Gaussian point
                        bgp = [0; -rho0*g]; % input body force
                        if norm(bgp) > 0
                            Fbeint = Fbeint + wtri*Ntri(6*(i-1)+2*gp-1:6*(i-1)+2*gp,:)'*bgp*latri(i)*h^2/4;
                        end
                    end
                end
                Ke = Keint; % element stiffness matrix
                Fbe = Fbeint; % element body force vector
            else % evaluate using exterior polyhedron
                Keext = zeros(8,8); % initialize element stiffness matrix for exterior domain of partial element
                Fbeext = zeros(8,1); % initialize element body force vector for exterior domain of partial element
                for i = 1:ntri % summation over triangles
                    for gp = 1:3 % summation over Gaussian points of triangles
                        Keext = Keext + wtri*Btri(9*(i-1)+3*gp-2:9*(i-1)+3*gp,:)'*D*Btri(9*(i-1)+3*gp-2:9*(i-1)+3*gp,:)*latri(i)*h^2/4;
                        %xgp = gcpolgp(1,3*(i-1)+gp); % x coordinate of Gaussian point
                        %ygp = gcpolgp(2,3*(i-1)+gp); % y coordinate of Gaussian point
                        bgp = [0; -rho0*g]; % input body force
                        if norm(bgp) > 0
                            Fbeext = Fbeext + wtri*Ntri(6*(i-1)+2*gp-1:6*(i-1)+2*gp,:)'*bgp*latri(i)*h^2/4;
                        end
                    end
                end
                Ke = Kefull-Keext; % element stiffness matrix
                Fbe = Fbefull-Fbeext; % element body force vector
            end
            Fte = zeros(8,1); % initialize element traction vector
            for i = 1:nseg % summation over segments
                for gp = 1:2 % summation over Gaussian points of segments
                    tgp = [0; 0]; % input traction
                    if norm(tgp) > 0
                        Fte = Fte + wseg*Nseg(4*(i-1)+2*gp-1:4*(i-1)+2*gp,:)'*tgp*llseg(i)*h/2;
                    end
                end
            end
        end
    end
end

function solid_area_percentage = find_solid_area(v1)
void_area = pi*v1(26:2:75)*v1(27:2:75)';
total_area = 2*2;
solid_area_percentage = (total_area - void_area)/total_area;
disp('solid area percentage is : ')
disp(solid_area_percentage)
end