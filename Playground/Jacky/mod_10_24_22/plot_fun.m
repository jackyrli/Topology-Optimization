function figure_handle = plot_fun(filename)
% Elastic moduli
lambda = 100e6; % elastic modulus
mu = 40e6; % elastic modulus
rho0 = 0e3; % density
g = 9.8; % g
displ = -1e-3; % displacement of the top surface

% Meshing
nex = 300; % number of elements in x direction
ney = 300; % number of elements in y direction
LX = 2; % length of mesh domain
LY = 2; % length of mesh domain
Xmin = 0; % minimum X coordinate
Ymin = 0; % minimum Y coordinate
Xmax = Xmin+LX; % maximum X coordinate
Ymax = Ymin+LY; % maximum Y coordinate
h = LX/nex; % mesh size

% Display settings
displayerror = 1; % whether display error message
displayface = 1; % select face color; 1: interp, 2: texturemap
displayedge = 0; % select edge color; 0: none, 1: white, 2: black
hidenodes = 0; % hide nodes near singularities

% Input file
vdata = load(filename);
v1 = zeros(25,1);
v2 = zeros(50,1);
v3 = zeros(25,1);

v1 = vdata(1:25);
v2 = vdata(26:75)+1e-15;
v3 = vdata(76:100);

gi = find(v1==1); % non-zero grid indices

% Elliptical boundaries
ngx = 5; % number of grids in x direction
ngy = 5; % number of grids in y direction
nN = length(gi); % number of Neumann ellipses
XN = Xmin+h*nex/ngx*(floor((gi-1)/ngy)+1/2); % x coordinates of centers of Neumann ellipses
YN = Ymin+h*ney/ngy*(mod(gi-1,ngy)+1/2); % y coordinates of centers of Neumann ellipses
AN = v2(2*gi-1); % semi-major axes of Neumann ellipses
BN = v2(2*gi); % semi-minor axes of Neumann ellipses
thetaN = v3(gi); % inclined angles of Neumann ellipses

% Basic parameters
nnx = nex+1; % number of nodes in x direction
nny = ney+1; % number of nodes in y direction
ne = nex*ney; % total number of elements
nn = (nex+1)*(ney+1); % total number of nodes
dof = 2*nn; % total degrees of freedom
lcgpe = 1/sqrt(3)*[-1 -1; 1 -1; -1 1; 1 1]'; % Gaussian points in local coordinates
we = 1; % Gaussian weight
epsilon = 1e-12; % epsilon for crossing points
maxit = 32; % number of maximum iterations
alpha = 1e7; % constant to normalize Dirichlet BC

if abs(LX/nex-LY/ney) > epsilon && displayerror == 1
    error('Error in meshing. Elements must be square-shaped.')
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% MESHING

% Connectivity array
connectivity = zeros(ne,4); % initialize connectivity array
for n = 1:4 % loop over element nodes
    for i = 1:nex
        for j = 1:ney
            connectivity(i+nex*(j-1),n) = mod(n-1,2)+i+(nex+1)*(floor((n-1)/2)+j-1); % insert node number
        end
    end
end

% Coordinates array
coordinates = zeros(2,nn); % initialize coordinates array
for i = 0:nex
    for j = 0:ney
        coordinates(:,1+i+(nex+1)*j) = [Xmin; Ymin] + h*[i;j]; % insert coordinate
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% LOCAL ARRAYS

% Interpolation arrays
Ne = zeros(8,8); % initialize [Ne] at Gaussian points
Be = zeros(12,8); % initialize [Be] at Gaussian points
for gp = 1:4 % loop over Gaussian points
    for n = 1:4 % loop over element nodes
        xigp = lcgpe(1,gp);
        etagp = lcgpe(2,gp);
        Ne(2*gp-1:2*gp,2*n-1:2*n) = 1/4*(1+(-1)^n*xigp)*(1+(-1)^floor((n+1)/2)*etagp)*eye(2);
        Be(3*gp-2,2*n-1) = 2/h*1/4*(-1)^n*(1+(-1)^floor((n+1)/2)*etagp);
        Be(3*gp-1,2*n) = 2/h*1/4*(-1)^floor((n+1)/2)*(1+(-1)^n*xigp);
        Be(3*gp,2*n-1) = 2/h*1/4*(-1)^floor((n+1)/2)*(1+(-1)^n*xigp);
        Be(3*gp,2*n) = 2/h*1/4*(-1)^n*(1+(-1)^floor((n+1)/2)*etagp);
    end
end

% Element elasticity matrix
D = [lambda+2*mu lambda 0; lambda lambda+2*mu 0; 0 0 mu]; % [D] (3x3)

% Element stiffness matrix
Kefull = zeros(8,8); % initialize [Ke] (8x8)
for gp = 1:4 % summation over Gaussian points
    Kefull = Kefull + we*Be(3*gp-2:3*gp,:)'*D*Be(3*gp-2:3*gp,:)*h^2/4;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% GLOBAL ARRAYS

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

disp('Assembling the global system...');
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
            xgp = gcgp(1,gp); % x coordinate of Gaussian point
            ygp = gcgp(2,gp); % y coordinate of Gaussian point
            bgp = [0; -rho0*g]; % input body force
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
                    xigp = (xgp-xc)*2/h;
                    etagp = (ygp-yc)*2/h;
                    Nb(2*gp-1:2*gp,2*n-1:2*n) = 1/4*(1+(-1)^n*xigp)*(1+(-1)^floor((n+1)/2)*etagp)*eye(2);
                end
                tgp = [0; 0]; % input traction
                if norm(tgp) > 0
                    Ftefull = Ftefull + we*Nb(2*gp-1:2*gp,:)'*tgp*h/2;
                end
            end
        end
        if xmax == Xmax % input full element Neumann boundary equation
            xgp = xmax; % x coordinate of Gaussian point
            Nb = zeros(4,8); % initialize [Nb] at Gaussian points
            for gp = 1:2 % loop over Gaussian points
                ygp = gcgp(2,2*gp-1); % y coordinate of Gaussian point
                for n = 1:4 % loop over element nodes
                    xigp = (xgp-xc)*2/h;
                    etagp = (ygp-yc)*2/h;
                    Nb(2*gp-1:2*gp,2*n-1:2*n) = 1/4*(1+(-1)^n*xigp)*(1+(-1)^floor((n+1)/2)*etagp)*eye(2);
                end
                tgp = [0; 0]; % input traction
                if norm(tgp) > 0
                    Ftefull = Ftefull + we*Nb(2*gp-1:2*gp,:)'*tgp*h/2;
                end
            end
        end
        if ymin == Ymin % input full element Neumann boundary equation
            ygp = ymin;
            Nb = zeros(4,8); % initialize [Nb] at Gaussian points
            for gp = 1:2 % loop over Gaussian points
                xgp = gcgp(1,gp); % x coordinate of Gaussian point
                for n = 1:4 % loop over element nodes
                    xigp = (xgp-xc)*2/h;
                    etagp = (ygp-yc)*2/h;
                    Nb(2*gp-1:2*gp,2*n-1:2*n) = 1/4*(1+(-1)^n*xigp)*(1+(-1)^floor((n+1)/2)*etagp)*eye(2);
                end
                tgp = [0; 0]; % input traction
                if norm(tgp) > 0
                    Ftefull = Ftefull + we*Nb(2*gp-1:2*gp,:)'*tgp*h/2;
                end
            end
        end
        if ymax == Ymax % input full element Neumann boundary equation
            ygp = ymax;
            Nb = zeros(4,8); % initialize [Nb] at Gaussian points
            for gp = 1:2 % loop over Gaussian points
                xgp = gcgp(1,gp); % x coordinate of Gaussian point
                for n = 1:4 % loop over element nodes
                    xigp = (xgp-xc)*2/h;
                    etagp = (ygp-yc)*2/h;
                    Nb(2*gp-1:2*gp,2*n-1:2*n) = 1/4*(1+(-1)^n*xigp)*(1+(-1)^floor((n+1)/2)*etagp)*eye(2);
                end
                tgp = [0; 0]; % input traction
                if norm(tgp) > 0
                    Ftefull = Ftefull + we*Nb(2*gp-1:2*gp,:)'*tgp*h/2;
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
            if abs(gcen(1,n)-Xmin) < epsilon || abs(gcen(1,n)-Xmax) < epsilon || abs(gcen(2,n)-Ymin) < epsilon || abs(gcen(2,n)-Ymax) < epsilon % input full element Dirichlet boundary equation
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
            [Ke,Fbe,Fte,typen] = NeumannBC(XN(Nn),YN(Nn),AN(Nn),BN(Nn),thetaN(Nn),h,e,en,gcen,xmin,ymin,xmax,ymax,nint,next,intn,extn,typen,D,rho0,g,Kefull,Fbefull,maxit,epsilon,displayerror); % Neumann partial element arrays
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
        if mod(percent,10) == 0
            disp([num2str(percent) '% assembly is completed.'])
        end
    end
    
end

K = sparse(spi,spj,spv,dof,dof); % form global stiffness matrix
K0 = K; % global stiffness matrix prior to imposing Dirichlet BC

toc
disp('Imposing Dirichlet boundary conditions...')
tic

% Impose Dirichlet BC on full elements
Dfulln = setdiff(Dfulln,0); % Dirichlet full element nodes
Dfulldof = zeros(2*length(Dfulln),1);
Dfullubar = zeros(2*length(Dfulln),1);
if ~isempty(Dfulln) % Dirichlet full element nodes are non-empty
    for i = 1:length(Dfulln) % loop over Dirichlet full element nodes
        n = Dfulln(i);
        %xn = coordinates(1,Dfulln(i));
        %yn = coordinates(2,Dfulln(i));
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

toc
disp('Prepare for solving...')
tic

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

%% SOLUTIONS AND PLOTS

toc
disp('Solving the global system...')
tic

% Solutions
uact = Kact\Fact; % solve for active nodes

toc
disp('Finished solving.')

% Process solutions
u = zeros(dof,1);
u(actdof) = uact;
u(sort([2*comn-1 2*comn])) = NaN;
r = zeros(dof,1);
r(actdof) = K0act*uact;
r(sort([2*comn-1 2*comn])) = NaN;

% Select plotting styles
switch displayface
    case 1
        surfface = 'interp';
    case 2
        surfface = 'texturemap';
    otherwise
        if displayerror == 1
            error('Invalid selection of face color.')
        end
end
switch displayedge
    case 0
        surfedge = 'none';
    case 1
        surfedge = 'w';
    case 2
        surfedge = 'k';
    otherwise
        if displayerror == 1
            error('Invalid selection of edge color.')
        end
end

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
    epsilon11(e) = (u1(en(2))+u1(en(4))-u1(en(1))-u1(en(3)))/(2*h);
    epsilon22(e) = (u2(en(3))+u2(en(4))-u2(en(1))-u2(en(2)))/(2*h);
    epsilon12(e) = (u1(en(3))+u1(en(4))-u1(en(1))-u1(en(2))+u2(en(2))+u2(en(4))-u2(en(1))-u2(en(3)))/(4*h);
    sigma11(e) = lambda*(epsilon11(e)+epsilon22(e)+epsilon33(e))+2*mu*epsilon11(e);
    sigma22(e) = lambda*(epsilon11(e)+epsilon22(e)+epsilon33(e))+2*mu*epsilon22(e);
    sigma33(e) = lambda*(epsilon11(e)+epsilon22(e)+epsilon33(e))+2*mu*epsilon33(e);
    sigma12(e) = 2*mu*epsilon12(e);
    sigma23(e) = 2*mu*epsilon23(e);
    sigma31(e) = 2*mu*epsilon31(e);
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
R2top = F2top-F2row/2
R2bottom = F2bottom-F2row/2
FEMweight=R2top+R2bottom;
% Weight
Vrect = LX*LY;
Vvoid = pi*dot(AN,BN);
Vmatl = Vrect-Vvoid;
weight = rho0*g*Vmatl;
% Stiffness
stiffnesspervolume = (abs(R2top)+abs(R2bottom))/2/abs(displ)/Vmatl
disp(stiffnesspervolume)
% von Mises stress
nactn = nex*ney-length(find(isnan(sigmavM)));
avgvM = sum(sigmavM,'omitnan')/nactn;
maxvM = max(sigmavM);
minvM = min(sigmavM);

% Display results
%disp(u)
%disp(typen)
%disp(avgvM)
%disp(maxvM)
%disp(minvM)
%disp(['Average von Mises stress = ',num2str(avgvM),' Pa'])
%writematrix(avgvM,'average_stress')

% Surface plots interpolated from nodal solutions
xn = reshape(coordinates(1,:),nnx,nny)';
yn = reshape(coordinates(2,:),nnx,nny)';
xc = repmat(linspace(Xmin+h/2,Xmax-h/2,nex),[ney 1]);
yc = repmat(linspace(Ymin+h/2,Ymax-h/2,ney)',[1 nex]);
u1n = reshape(u1,nnx,nny)';
u2n = reshape(u2,nnx,nny)';
r1n = reshape(r1,nnx,nny)';
r2n = reshape(r2,nnx,nny)';
sn = reshape(sn,nnx,nny)';
sigma11c = reshape(sigma11,nex,ney)';
sigma22c = reshape(sigma22,nex,ney)';
sigma12c = reshape(sigma12,nex,ney)';
sigmavMc = reshape(sigmavM,nex,ney)';
sigma11c([1:hidenodes ney-hidenodes+1:ney],[1:hidenodes nex-hidenodes+1:nex]) = NaN;
sigma22c([1:hidenodes ney-hidenodes+1:ney],[1:hidenodes nex-hidenodes+1:nex]) = NaN;
sigma12c([1:hidenodes ney-hidenodes+1:ney],[1:hidenodes nex-hidenodes+1:nex]) = NaN;
sigmavMc([1:hidenodes ney-hidenodes+1:ney],[1:hidenodes nex-hidenodes+1:nex]) = NaN;
figure_handle = figure(1)
subplot(2,3,1)
surf(xn,yn,u1n,'facecolor',surfface,'edgecolor',surfedge)
view(2)
axis equal
axis([Xmin Xmin+LX Ymin Ymin+LY])
title('u_1')
xlabel('x')
ylabel('y','rotation',0)
grid off
colorbar
subplot(2,3,2)
surf(xn,yn,u2n,'facecolor',surfface,'edgecolor',surfedge)
view(2)
axis equal
axis([Xmin Xmin+LX Ymin Ymin+LY])
title('u_2')
xlabel('x')
ylabel('y','rotation',0)
grid off
colorbar
subplot(2,3,3)
surf(xn,yn,sn,'facecolor',surfface,'edgecolor',surfedge)
view(2)
axis equal
axis([Xmin Xmin+LX Ymin Ymin+LY])
title('$$\mathbf{\sqrt{u_1^2+u_2^2}}$$','interpreter','latex')
xlabel('x')
ylabel('y','rotation',0)
grid off
colorbar
subplot(2,3,4)
%surf(xc,yc,sigma11c,'facecolor',surfface,'edgecolor',surfedge)
surf(xn,yn,r1n,'facecolor',surfface,'edgecolor',surfedge)
view(2)
axis equal
axis([Xmin Xmin+LX Ymin Ymin+LY])
%title('\sigma_{11}')
title('R_1')
xlabel('x')
ylabel('y','rotation',0)
grid off
colorbar
subplot(2,3,5)
%surf(xc,yc,sigma22c,'facecolor',surfface,'edgecolor',surfedge)
surf(xn,yn,r2n,'facecolor',surfface,'edgecolor',surfedge)
view(2)
axis equal
axis([Xmin Xmin+LX Ymin Ymin+LY])
%title('\sigma_{22}')
title('R_2')
xlabel('x')
ylabel('y','rotation',0)
grid off
colorbar
subplot(2,3,6)
surf(xc,yc,sigmavMc,'facecolor',surfface,'edgecolor',surfedge)
view(2)
axis equal
axis([Xmin Xmin+LX Ymin Ymin+LY])
title('von Mises Stress')
xlabel('x')
ylabel('y','rotation',0)
grid off
colorbar
sgtitle('Cut-Cell Finite Element Method for 2D Isotropic Linear Elastic Material with Traction-Free Ellipses')
set(gcf,'position',get(0,'screensize'));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% NEUMANN BC FOR PARTIAL ELEMETNS

function [Ke,Fbe,Fte,typen] = NeumannBC(XN,YN,AN,BN,thetaN,h,e,en,gcen,xmin,ymin,xmax,ymax,nint,next,intn,extn,typen,D,rho0,g,Kefull,Fbefull,maxit,epsilon,displayerror)

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