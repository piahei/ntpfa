function coSet = getCollactionSetOPT(G, rock)
if nargin == 1
    rock = makeRock(G, 1, 1);
end

n_hf = size(G.cells.faces, 1);
n_f = G.faces.num;
dim = G.griddim;
numFace = diff(G.cells.facePos);
maxNumFace = max(numFace); %maximum number of faces in cells
isBF = false(G.faces.num, 1);
isBF(boundaryFaces(G)) = true;
n_bf = sum(isBF);
n_if = n_f-n_bf;
isIF = (1:1:n_f)';
isIF = isIF(~isBF);

faceSign = zeros(n_hf, 1);
N = G.faces.neighbors;
% This builds a collocation set for the nonlinear two-point flux
% approximation. Generally, we construct the necessary coefficients
% and vectors using optimization as a preprocessing step.
%
% All of these sets are constructed twice for internal faces

% Vectors from d.o.f. to interpolation points
T_all = cell(maxNumFace, 2);
[T_all{:}] = deal(zeros(n_f, dim));

% Coefficient matrices for each set of vectors
[coeff_ij, coeff_ji] = deal(zeros(n_if, maxNumFace));
coeff_bf = zeros(n_bf,maxNumFace);

% Active indicators for each set of vectors/coefficients. A point is
% said to be active if it is defined in the cell where we want a
% degree of freedom for a two point method (plus the face between two
% cells).
[act_ij, act_ji] = deal(false(n_if, maxNumFace));
act_cell_bf = false(n_bf,maxNumFace);


% Cell neighbors for all cells needed when computing the flux.

[cn_1, cn_2]     = deal(ones(n_if, maxNumFace));
cn_bf = ones(n_bf,maxNumFace);

% L is the face normal vectors, scaled by cell permeability and face
% area. One for each half face.
L = zeros(n_hf, dim);
L_face = zeros(n_f, dim);

% cellNo = rldecode(1 : G.cells.num, ...
%     diff(G.cells.facePos), 2) .';
% faceNo = G.cells.faces(:, 1);

%Find variables for each face
[hap,~,~,O_i,O_j] = findVariablesNTPFA(G,rock);

% Preallocate matrices
L1_ij = zeros(n_if,maxNumFace);
L2_ij = zeros(n_if,maxNumFace);
L1_ji = zeros(n_if,maxNumFace);
L2_ji = zeros(n_if,maxNumFace);
A_ij = zeros(n_if,maxNumFace);
A_ji = zeros(n_if,maxNumFace);

L1_bf = zeros(n_bf,maxNumFace);
L2_bf = zeros(n_bf,maxNumFace);
A_bf = zeros(n_bf,maxNumFace);

bf = 0;
iF = 0;
for i = 1:n_if
    f = isIF(i);
%     if isBF(f)
%         bf=bf+1;
%         c = N(f,:); c(c==0)=[];
%         nF = numFace(c);
%         sgn = 1 - 2*(N(f, 2) == c);
%         faceSign(i) = sgn;
%         normal = G.faces.normals(f, :);
%         K = getK(rock, c, dim);
%         l = normal*sgn*K;
%         L(i, :) = l;
%         i = i+1;
%         L_face(f, :) = normal*K;
%         cpt = G.cells.centroids(c,:);
%         selfFaces = gridCellFaces(G,c);
%         intSelfFaces = selfFaces(~isBF(selfFaces));
%         nIF = numel(intSelfFaces);
%         
%         [T, ~, coeff, ~, successCell] = collocateSetOPT(G,cpt,intSelfFaces,hap(intSelfFaces,:), l);
%         
%         assert(successCell)
%         
%         for j = 1:nIF
%             T_all{j, 1}(f, :) = T(j, :);  %brukes ikke
%         end
%         
%         coeff_bf(bf,1:nIF) = coeff;
%         n = N(selfFaces,:); n=sum(n,2)-c; %neighbors to cell, including "ouside" cells
%         act_cell_bf(bf,~isBF(selfFaces)) = true; %active cell neighbors
%         n(n==0)=1; %set outside cell to 1 to avoid nonlogical indexing
%         cn_bf(bf,1:nF) = n;
%         A_bf(bf,1:nIF) = coeff.*var.O_j(intSelfFaces);
%         L1_bf(bf,1:nIF) = A_bf(bf,1:nIF)'.*var.O_i(intSelfFaces);
%         L2_bf(bf,1:nIF) = A_bf(bf,1:nIF)'.*var.O_j(intSelfFaces);
%         
%     else
        
        iF = iF+1;
        c1 = N(f,1);
        nF_c1 = numFace(c1); 
        c2 = N(f,2);
        nF_c2 = numFace(c2); 
        
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
        sgn = 1 - 2*(G.faces.neighbors(f, 2) == c1); %unødvendig, kan settes til +1
        faceSign(i) = sgn; %brukes ikke..
        faceSign(i+1) = -sgn;
        normal = G.faces.normals(f, :);
        K = getK(rock, c1, dim);
        l = normal*sgn*K;
        L(i, :) = l;
        L(i+1,:) = -l;
        L_face(f, :) = normal*K;
        
        % Collocation accross face direction ij
        cpt = G.cells.centroids(c1, :);  %centerpoint of cell 1
        
        [T, ~, coeff, ~, successCell] = collocateSetOPT(cpt,selfFaces1,hap(selfFaces1,:), l);

        assert(successCell)
        
        A_ij(iF,1:nF_c1) = coeff.*O_j(selfFaces1);  
        L1_ij(iF,1:nF_c1) = A_ij(iF,1:nF_c1)'.*O_i(selfFaces1);
        L2_ij(iF,1:nF_c1) = A_ij(iF,1:nF_c1)'.*O_j(selfFaces1);
         
        
        for j = 1:nF_c1
            T_all{j, 1}(f, :) = T(j, :);
        end
        
        coeff_ij(iF,1:nF_c1) = coeff;
       
        % Collocation accross face direction ji
        
        cpt = G.cells.centroids(c2,:);
      
        [T, ~, coeff, ~, successFace] = collocateSetOPT(cpt, selfFaces2, hap(selfFaces2,:), -l);
        
        assert(successFace)
        
        A_ji(iF,1:nF_c2) = coeff.*O_j(selfFaces2);
        L1_ji(iF,1:nF_c2) = A_ji(iF,1:nF_c2)'.*O_i(selfFaces2);
        L2_ji(iF,1:nF_c2) = A_ji(iF,1:nF_c2)'.*O_j(selfFaces2);
        
        
        for j = 1:nF_c2
            T_all{j, 2}(f, :) = T(j, :);
        end
        
        coeff_ji(iF,1:nF_c2) = coeff;
      
%     end
end
bndryVars = struct('L1_bf',L1_bf,'L2_bf',L2_bf,'A_bf',A_bf);
intVars = struct('L1_ij',L1_ij,'L2_ij',L2_ij,'L1_ji',L1_ji,'L2_ji',L2_ji,...
    'A_ij',A_ij,'A_ji',A_ji);
vars = struct('internalVariables',intVars, 'boundaryVariables',bndryVars);
coSet.variables = vars;
%     coSet.L_ij = {L1_ij,L2_ij};
%     coSet.L_ji = {L1_ji,L2_ji};
%     coSet.L_bf = L_bf;
%     coSet.A = {A_ij,A_ji,A_bf}

coSet.T = T_all;
%     coSet.T_norm = cellfun(@(x) sqrt(sum(x.^2, 2)), T_all, 'UniformOutput', false);
C = struct('coeff_ij',coeff_ij,'coeff_ji',coeff_ji);
coSet.C = C;
coSet.cellNeighbors = {cn_1, cn_2, cn_bf}; 
coSet.l = L;
coSet.active = struct('act_ij',act_ij,'act_ji',act_ji,...
    'act_cell_bf',act_cell_bf);

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

function [T, types, coefficients, globIx, ok] = collocateSetOPT(cpt,faces,hap, l)
% [pts, typ, glob] = getCandidatePointsOPT(G,faces);
warning('off','MATLAB:singularMatrix');
lb = zeros(length(hap)+1,1);% lb(1:end-1)=-1;
[T, ix, coefficients, ok] = OPTgetSetForHF(cpt, l, hap,lb,0);

globIx = faces(ix);
types = ones(4,1)*2;
end

function [t, ixSet, coeff, done] = OPTgetSetForHF(center,l,hap,lb,it)
% lb is lower boundary for optimization. Input is given as 0,
% and is lowered if valid positive decomposition does not exist.

% quit if valid decomposition is not found after 40 iterations 
if it == 40
    done = false;
    error('no decomposition found after 40 iterations')
end

ixSet = (1:1:length(hap));
%dim = size(pts,2);
N = size(hap,1);  %number of faces in current cell
coeff = zeros(N+1,1);

t = bsxfun(@minus,hap,center);
t_v = sqrt(sum(t.^2,2));   %norm of each row in t
t_u = bsxfun(@rdivide, t, t_v)';  %normalized t
T = (t_u(:,ixSet))'; %for debugging

l_v = norm(l,2);  %norm of facenormal
l_u = l./l_v;     %normalized facenormal

s = bsxfun(@rdivide,t_v,l_v);  %sigma
k = N/min(s)+10000;            %kappa >> N/min(s)
if isnan(k) || isinf(k)
    k = 10000;
end

%vector for objective function. Last element is gamma
c = ones(N+1,1);
c(end) = k;

funk = @(x) c'*x;
%add zeros to end of t-vector to allow for gamma coefficient. 
t_u(:,N+1)=0;

% Set tolerance for optimization, nonlinear constraints are given in @con
% Set lb to 0 to begin with, lower if valid positive decompotition does not exist

lb(1:end-1)=-1000;

option = optimset('Display','off', 'TolCon',1e-20);%,'Algorithm','sqp');
coeff = fmincon(funk,coeff,[],[], t_u, l_u',lb,[],@(coeff)con(coeff,s),option);
g = coeff(end);
coeff(end) = [];
done = true;
% check sign on decomposition. Set sign(0) to 1 to give slack.
% sgn_l_u = sign(l_u)'; 
% sgn_decomp = sign(T'*coeff);
% if any(sgn_l_u == 0) 
%     sgn_l_u(sgn_l_u == 0) = 1;
% end
% if any(sgn_decomp == 0)
%     sgn_decomp(sgn_decomp==0)=1;
% elseif any(abs(T'*coeff)<1e-5)
%     sgn_decomp(abs(T'*coeff)<1e-5) = 1; %round off 0
% end
% 
% % check criteria. Adjust lb if criteria is not met. 
% if all(sgn_decomp == sgn_l_u) && all(coeff >= -g) && sum(coeff./s)>=0 
%     done = true;
% else
%     it = it+1;
%     if it==1
%         lb(1:end-1) = -0.1; %lower lb
%     else
%         lb = lb*1.5;  %lower lb
%     end
%     [T, ixSet, coeff, done] = OPTgetSetForHF(center,l,hap,lb,it);
% end
coeff = coeff./s;

end

function [c,ceq] = con(x,s)
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
        centroid1 = G.cells.centroids(c1,:); %centroid of neighboring cell 1
        centroid2 = G.cells.centroids(c2,:); %centroid of neighboring cell 2
        dist_i = dist(fn1,fn2,centroid1);    % shortest distance from cell centroid 1 to face
        dist_j = dist(fn1,fn2,centroid2);    % shortest distance from cell centroid 2 to face
        beta_ij(f) = (G.faces.normals(f,:)*K1*G.faces.normals(f,:)')/dist_i;
        beta_ji(f) = (G.faces.normals(f,:)*K2*G.faces.normals(f,:)')/dist_j;
        
        if ~isfinite(beta_ij(f)) %why inf?! SEAMOUNT TRIANGLE
            beta_ij(f) = 0;
        elseif ~isfinite(beta_ji(f))
            beta_ji(f) = 0;
        end
        hap(f,:) = (beta_ij(f)+beta_ji(f))\...
            (beta_ij(f)*centroid1'+beta_ji(f)*centroid2'+(K1-K2)*G.faces.normals(f,:)');
        if any(~isfinite(hap(f,:)))
            hap(f,:) = G.faces.centroids(f,:);
        end
    end %if
    O_i(f) = beta_ij(f)/(beta_ij(f)+beta_ji(f));  %omega
    O_j(f) = beta_ji(f)/(beta_ij(f)+beta_ji(f));
end %for

end %function

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
        centroid1 = G.cells.centroids(c1,:); %centroid of neighboring cell 1
        centroid2 = G.cells.centroids(c2,:); %centroid of neighboring cell 2
        dist_i = dist(fCntr,centroid1,normal);    % shortest distance from cell centroid 1 to face
        dist_j = dist(fCntr,centroid2,normal);    % shortest distance from cell centroid 2 to face
        beta_ij(f) = (fCntr*K1*fCntr')/dist_i;
        beta_ji(f) = (fCntr*K2*fCntr')/dist_j;
        
        if ~isfinite(beta_ij(f)) %why inf?! SEAMOUNT TRIANGLE
            beta_ij(f) = 0;
        elseif ~isfinite(beta_ji(f))
            beta_ji(f) = 0;
        end
        hap(f,:) = (beta_ij(f)+beta_ji(f))\...
            (beta_ij(f)*centroid1'+beta_ji(f)*centroid2'+(K1-K2)*G.faces.normals(f,:)');
        if any(~isfinite(hap(f,:)))
            hap(f,:) = fCntr;
        end
    end %if
    O_i(f) = beta_ij(f)/(beta_ij(f)+beta_ji(f));  %omega
    O_j(f) = beta_ji(f)/(beta_ij(f)+beta_ji(f));
end 

end 
