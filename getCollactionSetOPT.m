function coSet = getCollactionSetOPT(G, rock)
% This builds a collocation set for the nonlinear two-point flux
% approximation method using optimization to generate the necessary
% coefficients and vectors. Used in PressureOilWaterModelNTPFAopt.
% 
% SYNOPSIS:
%   soSet = getCollactionSetOPT(G, rock)
%
% PARAMETERS:
%   G    - Grid structure as described by grid_structure.
%
%   rock - Rock data structure with valid field `perm`.  The permeability
%          is assumed to be in measured in units of metres squared (m^2).
%          Use function `darcy` to convert from darcies to m^2, e.g::
%
%                 perm = convertFrom(perm, milli*darcy)
%
%          if the permeability is provided in units of millidarcies.
%
%          The field rock.perm may have ONE column for a scalar
%          permeability in each cell, TWO/THREE columns for a diagonal
%          permeability in each cell (in 2/3 D) and THREE/SIX columns for a
%          symmetric full tensor permeability.  In the latter case, each
%          cell gets the permeability tensor::
%
%                 K_i = [ k1  k2 ]      in two space dimensions
%                       [ k2  k3 ]
%
%                 K_i = [ k1  k2  k3 ]  in three space dimensions
%                       [ k2  k4  k5 ]
%                       [ k3  k5  k6 ]
% RETURNS:
%  coSet - collocation set for the NTPFA method using optimization for
%          decomposition
%
% SEE ALSO:
%   `computeGeometry`, `makeRock`, 'PressureOilWaterModelNTPFAopt'.

%{
Written by Pia-Kristina Heigrestad, 2018 
%}

if nargin == 1
    rock = makeRock(G, 1, 1);
end

n_f = G.faces.num;
dim = G.griddim;
numFace = diff(G.cells.facePos);
maxNumFace = max(numFace); %maximum number of faces in cells
isIF = true(G.faces.num, 1);
isIF(boundaryFaces(G)) = false;
n_if = sum(isIF);
faces = 1:1:n_f;
isIF = faces(isIF);

N = G.faces.neighbors;


% Vectors from d.o.f. to interpolation points
T_all = cell(maxNumFace, 2);
[T_all{:}] = deal(zeros(n_f, dim));

% Coefficient matrices. Each row corresponds to an internal face.
[coeff_ij, coeff_ji] = deal(zeros(n_if, maxNumFace));

% Active indicators for cell neighbors for each set of vectors/coefficients.
[act_ij, act_ji] = deal(false(n_if, maxNumFace));

% Cell neighbors for all cells needed when computing the flux.

[cn_1, cn_2]     = deal(ones(n_if, maxNumFace));

%Find variables for each face
[hap,~,~,O_i,O_j] = findVariablesNTPFA(G,rock);

% Preallocate matrices for coeff*omega
A_ij = zeros(n_if,maxNumFace);
A_ji = zeros(n_if,maxNumFace);

% All variables are found in a for-loop with length equal to 
% number of internal faces. Two sets for each face.
for iF = 1:n_if
    f = isIF(iF);
    c1 = N(f,1);
    cpt1 = G.cells.centroids(c1, :);  %centroid of cell 1
    nF_c1 = numFace(c1);  %number of faces in cell 1
    c2 = N(f,2);
    cpt2 = G.cells.centroids(c2,:);  %centroid of cell 2
    nF_c2 = numFace(c2);  %number of faces in cell 2
    
    % Self faces are faces belonging to current cell
    selfFaces1 = gridCellFaces(G, c1);
    n1 = N(selfFaces1,:); n1=sum(n1,2)-c1; %neighbors to cell, including "outside" cells
    n1(n1==0)=c1; %set 'outside pressure' equal to current cell for continuity
    cn_1(iF,1:nF_c1) = n1;
    act_ij(iF,selfFaces1~=f) = true;
    selfFaces2 = gridCellFaces(G, c2);
    n2 = N(selfFaces2,:); n2=sum(n2,2)-c2;
    n2(n2==0)=c2;
    cn_2(iF, 1:nF_c2) = n2;
    act_ji(iF,selfFaces2~=f) = true;
    
    % Compute permeability scaled normal vector
    normal = G.faces.normals(f, :);
    K = getK(rock, c1, dim);
    l = normal*K;    
    
    % Collocation by optimization accross face direction ij
    
    [T, ~, coeff, successCell] = collocateSetOPT(cpt1,hap(selfFaces1,:), l);
    
    assert(successCell)
    
    A_ij(iF,1:nF_c1) = coeff.*O_j(selfFaces1);
    
    for j = 1:nF_c1
        T_all{j, 1}(f, :) = T(j, :);
    end
    
    coeff_ij(iF,1:nF_c1) = coeff;
    
    % Collocation by optimization accross face direction j. Change sign of
    % normal vector. 
    
    [T, ~, coeff, successFace] = collocateSetOPT(cpt2, hap(selfFaces2,:), -l);
    
    assert(successFace)
    
    A_ji(iF,1:nF_c2) = coeff.*O_i(selfFaces2);
    
    for j = 1:nF_c2
        T_all{j, 2}(f, :) = T(j, :);
    end
    
    coeff_ji(iF,1:nF_c2) = coeff;
    
end

% Store all variables.

coSet.variables = struct('A_ij',A_ij,'A_ji',A_ji);
coSet.T = T_all;
coSet.C = struct('coeff_ij',coeff_ij,'coeff_ji',coeff_ji);
coSet.cellNeighbors = {cn_1, cn_2};
coSet.active = struct('act_ij',act_ij,'act_ji',act_ji);

end

function K = getK(rock, cell, dim)
k = rock.perm(cell, :);
switch numel(k)
    case 1
        % Scalar perm
        K = k;
    case dim
        % Diagonal tensor
        K = diag(k);
    case 3*(dim - 1)
        % Full symmetric tensor
        if dim == 2
            K = [k(1), k(2); ...
                k(2), k(3)];
        else
            K = [k(1), k(2), k(3); ...
                k(2), k(4), k(5); ...
                k(3), k(5), k(6)];
        end
    otherwise
        error('What sorcery is this?!');
end
end

function [t, set, coeff, done] = collocateSetOPT(center,hap,l)
warning('off','MATLAB:singularMatrix');

set = (1:1:length(hap));
N = size(hap,1);  %number of faces in current cell
coeff = zeros(N+1,1);

t = bsxfun(@minus,hap,center);

t_v = sqrt(sum(t.^2,2));   %norm of each row in t
t_u = bsxfun(@rdivide, t, t_v)';  %normalized t

T = (t_u(:,set))'; %for debugging

l_v = norm(l,2);  %norm of facenormal
l_u = l./l_v;     %normalized facenormal

s = bsxfun(@rdivide,t_v,l_v);  %sigma
k = N/min(s)+10000;            %kappa >> N/min(s)
if isnan(k) || isinf(k)  %in case min(s)=0
    k = 10000;
end

%vector for objective function. Last element is gamma.
c = ones(N+1,1);
c(end) = k;

funk = @(x) c'*x;
%add zeros to end of t-vector to allow for gamma coefficient.
t_u(:,N+1)=0;

% Set tolerance for optimization, nonlinear constraints are given in @con
lb = zeros(length(hap)+1,1);
lb(1:end-1)=-1000;  %lower boundary
option = optimset('Display','off', 'TolCon',1e-20);

coeff = fmincon(funk,coeff,[],[], t_u, l_u',lb,[],@(coeff)con(coeff,s),option);
g = coeff(end);
coeff(end) = [];

% % check criteria. Display warning if not all constraints are met.
if all(coeff >= -g) && sum(coeff./s)>=0
    done = true;
else
    done = true;
    disp('Warning: fmincon did not succeed on meeting all constraints. Model will still be built.')
end

coeff = coeff./s;

end

function [c,ceq] = con(x,s)
%function for nonliear constraint
A = zeros(length(x));
A(1,1:end-1) = s.^(-1);
A(2:end,end) = 1;
A(2:end,1:end-1)=eye(length(x)-1);
c = -A*x;
ceq = [];
end

function [hap,beta_ij,beta_ji,O_i,O_j] = findVariablesNTPFA(G,rock)
dim = G.griddim;
if dim == 2
    [hap,beta_ij,beta_ji,O_i,O_j] = vars2D(G,rock);
elseif dim == 3
    [hap,beta_ij,beta_ji,O_i,O_j] = vars3D(G,rock);
end
end

function [hap,beta_ij,beta_ji,O_i,O_j] = vars2D(G,rock)
dim = 2;
noFaces = G.faces.num;
bFace = boundaryFaces(G);
N = G.faces.neighbors;

nodePos = G.faces.nodePos;
nodes = G.faces.nodes;

% Preallocate variables
hap = zeros(noFaces,dim); %harmonic average point
hap(bFace,:) = G.faces.centroids(bFace,:); %hap for bf is equal to centroid

beta_ij = zeros(noFaces,1);
beta_ji = zeros(noFaces,1);
O_i = zeros(noFaces,1); %omega_ij
O_j = zeros(noFaces,1); %omega_ji

% function for distance from cell centroid to face
dist = @(p1,p2,x) abs((p2(2)-p1(2))*x(1)-(p2(1)-p1(1))*x(2)+p2(1)*p1(2)-p2(2)*p1(1))/...
    sqrt((p2(2)-p1(2))^2+(p2(1)-p1(1))^2);

for f = 1:noFaces  %over all faces
    
    c1 = N(f,1); %cell neighbor 1
    c2 = N(f,2); %cell neighbor 2
    
    node = nodes(nodePos(f) : nodePos(f+1)-1, :);
    
    fn1 = G.nodes.coords(node(1),:);  %coordinates to node 1
    fn2 = G.nodes.coords(node(2),:);  %coordinates to node 2
    
    % for boundary faces, permeability is mirrored outside boundary
    if c1 == 0
        K2 = getK(rock,c2,dim);
        %  K1 = K2;
        centroid2 = G.cells.centroids(c2,:);
        dist_j = dist(fn1,fn2,centroid2);
        beta_ji(f) = (G.faces.normals(f,:)*K2*G.faces.normals(f,:)')/dist_j;
        if ~isfinite(beta_ji(f))
            beta_ji(f) = 0;
        end
        beta_ij(f) = beta_ji(f);
        
    elseif c2 == 0
        K1 = getK(rock,c1,dim);
        %  K2 = K1;
        centroid1 = G.cells.centroids(c1,:);
        dist_i = dist(fn1,fn2,centroid1);
        beta_ij(f) = (G.faces.normals(f,:)*K1*G.faces.normals(f,:)')/dist_i;
        if ~isfinite(beta_ij(f))
            beta_ij(f) = 0;
        end
        beta_ji(f) = beta_ij(f);
    else
        K1 = getK(rock,c1,dim);
        K2 = getK(rock,c2,dim);
        centroid1 = G.cells.centroids(c1,:); 
        centroid2 = G.cells.centroids(c2,:);
        dist_i = dist(fn1,fn2,centroid1);    
        dist_j = dist(fn1,fn2,centroid2);   
        beta_ij(f) = (G.faces.normals(f,:)*K1*G.faces.normals(f,:)')/dist_i;
        beta_ji(f) = (G.faces.normals(f,:)*K2*G.faces.normals(f,:)')/dist_j;
        
        %debugging for PEBI grid
        if ~isfinite(beta_ij(f)) 
            beta_ij(f) = 0;
        elseif ~isfinite(beta_ji(f))
            beta_ji(f) = 0;
        end
      
        hap(f,:) = (beta_ij(f)+beta_ji(f))\...
            (beta_ij(f)*centroid1'+beta_ji(f)*centroid2'+(K1-K2)*G.faces.normals(f,:)');
        if any(~isfinite(hap(f,:)))
            hap(f,:) = G.faces.centroids(f,:);
        end
    end 
    O_i(f) = beta_ij(f)/(beta_ij(f)+beta_ji(f)); 
    O_j(f) = beta_ji(f)/(beta_ij(f)+beta_ji(f));
end 
end

function [hap,beta_ij,beta_ji,O_i,O_j] = vars3D(G,rock)
dim = 3;
noFaces = G.faces.num;
bFace = boundaryFaces(G);
N = G.faces.neighbors;

% Preallocate variables
hap = zeros(noFaces,dim);
hap(bFace,:) = G.faces.centroids(bFace,:); %hap for bf is equal to centroid

beta_ij = zeros(noFaces,1); %beta for internal faces
beta_ji = zeros(noFaces,1);
O_i = zeros(noFaces,1); %omega_ij
O_j = zeros(noFaces,1); %omega_ji

% function for distance from cell centroid to face
dist = @(x0,x1,n) abs(n(1)*(x1(1)-x0(1))+n(2)*(x1(2)-x0(2))+n(3)*(x1(3)-x0(3)))/norm(n);

for f = 1:noFaces  %over all faces
    
    c1 = N(f,1); %cell neighbor 1
    c2 = N(f,2); %cell neighbor 2
    normal = G.faces.normals(f,:);
    fCntr = G.faces.centroids(f,:);
    
    % for boundary faces, permeability is mirrored outside boundary
    if c1 == 0
        K2 = getK(rock,c2,dim);
        %  K1 = K2;
        centroid2 = G.cells.centroids(c2,:);
        dist_j = dist(fCntr,centroid2,normal);
        beta_ji(f) = (normal*K2*normal')/dist_j;
        if ~isfinite(beta_ji(f))
            beta_ji(f) = 0;
        end
        beta_ij(f) = beta_ji(f);
        
    elseif c2 == 0
        K1 = getK(rock,c1,dim);
        %  K2 = K1;
        centroid1 = G.cells.centroids(c1,:);
        dist_i = dist(fCntr,centroid1,normal);
        beta_ij(f) = (fCntr*K1*fCntr')/dist_i;
        if ~isfinite(beta_ij(f))
            beta_ij(f) = 0;
        end
        beta_ji(f) = beta_ij(f);
    else
        K1 = getK(rock,c1,dim);
        K2 = getK(rock,c2,dim);
        centroid1 = G.cells.centroids(c1,:);
        centroid2 = G.cells.centroids(c2,:); 
        dist_i = dist(fCntr,centroid1,normal);    
        dist_j = dist(fCntr,centroid2,normal);    
        beta_ij(f) = (fCntr*K1*fCntr')/dist_i;
        beta_ji(f) = (fCntr*K2*fCntr')/dist_j;
        
        % debugging for PEBI grid
        if ~isfinite(beta_ij(f))
            beta_ij(f) = 0;
        elseif ~isfinite(beta_ji(f))
            beta_ji(f) = 0;
        end
        hap(f,:) = (beta_ij(f)+beta_ji(f))\...
            (beta_ij(f)*centroid1'+beta_ji(f)*centroid2'+(K1-K2)*G.faces.normals(f,:)');
        if any(~isfinite(hap(f,:)))
            hap(f,:) = fCntr;
        end
    end 
    O_i(f) = beta_ij(f)/(beta_ij(f)+beta_ji(f));  %omega
    O_j(f) = beta_ji(f)/(beta_ij(f)+beta_ji(f));
end 
end 
