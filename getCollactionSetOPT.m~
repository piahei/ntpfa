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
    
    %jumpFlag = getJumpFlag(G, rock);
    faceSign = zeros(n_hf, 1);
    N = G.faces.neighbors;
    % This builds a collocation set for the nonlinear two-point flux
    % approximation. Generally, we construct the necessary coefficients 
    % and vectors using optimization as a preprocessing step. 
    % 
    % Most of these sets are constructed twice (internal faces): 
    % Once across the face ij and once across ji 
    
    % Vectors from d.o.f. to interpolation points
    T_all = cell(maxNumFace, 2);  
    [T_all{:}] = deal(zeros(n_f, dim));
    
    % Coefficient matrices for each set of vectors
    [coeff_ij, coeff_ji]       = deal(zeros(n_if, maxNumFace));
    coeff_bf = zeros(n_bf,maxNumFace);
    % Active indicators for each set of vectors/coefficients. A point is
    % said to be active if it is defined in the cell where we want a
    % degree of freedom for a two point method (plus the face between two
    % cells).
    [act_cell_ij, act_cell_ji, act_cell_nghbr_i, act_cell_nghbr_j] = deal(false(n_if, maxNumFace));
    act_cell_bf = false(n_bf,maxNumFace);
    
   
    % Global indices for the pressure values we need to sample for each
    % point. Cell neighbors for cells needed when computing the flux.
    
    [cn_1, cn_2]     = deal(zeros(n_if, maxNumFace));
    cn_bf = zeros(n_bf,maxNumFace);
    % Type array indicating what the indices actually point to:
    %         1: Cell
    %         2: Face
    %         3: Edge ??
    %[type_ij, type_ji] = deal(ones(n_f, numFace)*2);  Removed this
    
    % L is the face normal vectors, scaled by cell permeability and face
    % area. One for each half face.
    L = zeros(n_hf, dim);
    L_face = zeros(n_f, dim);
    
    cellNo = rldecode(1 : G.cells.num, ...
                     diff(G.cells.facePos), 2) .';  
    faceNo = G.cells.faces(:, 1);
   
    
    %Find variables for each face
    [hap,beta_ij,beta_ji,O_i,O_j] = findVariablesNTPFA(G,rock); 
    var = struct('beta_ij',beta_ij,'beta_ji',beta_ji,'O_i',O_i,'O_j',O_j);
    %coSet.variables = var;
    
    L1_ij = zeros(n_if,maxNumFace);
    L2_ij = zeros(n_if,maxNumFace);
    L1_ji = zeros(n_if,maxNumFace);
    L2_ji = zeros(n_if,maxNumFace);
    A_ij = zeros(n_if,maxNumFace);
    A_ji = zeros(n_if,maxNumFace);
    
    L1_bf = zeros(n_bf,maxNumFace);
    L2_bf = zeros(n_bf,maxNumFace);
    A_bf = zeros(n_bf,maxNumFace);
    
    i = 1;  %index for half face
    bf = 0;
    iF = 0;
    for f = 1:n_f 
        if isBF(f)
            bf=bf+1;
            c = N(f,find(N(f,:)~=0)); 
            nF = numFace(c);
            sgn = 1 - 2*(N(f, 2) == c);
            faceSign(i) = sgn;
            normal = G.faces.normals(f, :);
            K = getK(rock, c, dim);
            l = normal*sgn*K;
            L(i, :) = l;
            i = i+1;
            L_face(f, :) = normal*K;
            cpt = G.cells.centroids(c,:);
            selfFaces = gridCellFaces(G,c);
            
            [T, ~, coeff, ~, successCell] = collocateSetOPT(G,cpt,selfFaces,hap(selfFaces,:), l);
            
            assert(successCell)
            
            for j = 1:nF  
                T_all{j, 1}(f, :) = T(j, :);
            end
            
            coeff_bf(bf,1:nF) = coeff;
            n = N(selfFaces,:); n(n==c)=[]; %neighbors to cell, including "ouside" cells
            act_cell_bf(bf,find(n~=0)) = true; %active cell neighbors
            n(n==0)=1; %set outside cell to 1 to avoid nonlogical indexing
            cn_bf(bf,1:nF) = n;  
            A_bf(bf,1:nF) = coeff.*var.O_j(selfFaces);
            L1_bf(bf,1:nF) = A_bf(bf,1:nF)'.*var.O_i(selfFaces);
            L2_bf(bf,1:nF) = A_bf(bf,1:nF)'.*var.O_j(selfFaces);
            
        else
            
            iF = iF+1;
            c1 = N(f,1);  
            nF_c1 = numFace(c1);
            c2 = N(f,2);  
            nF_c2 = numFace(c2);
            % Self faces are faces belonging to current cell
            selfFaces1 = gridCellFaces(G, c1);
            n_i = find(selfFaces1==f);
            act_cell_nghbr_i(iF,n_i) = true;
            selfFaces2 = gridCellFaces(G, c2);
            n_j = find(selfFaces2==f);
            act_cell_nghbr_j(iF,n_j) = true;
            %adjFaces = adjacentFacesForFace(G, f, dim - 1);
            
            % Compute permeability scaled normal vector
            sgn = 1 - 2*(G.faces.neighbors(f, 2) == c1); %unødvendig, kan settes til +1
            faceSign(i) = sgn;
            faceSign(i+1) = -sgn;
            normal = G.faces.normals(f, :);
            K = getK(rock, c1, dim);
            l = normal*sgn*K;
            L(i, :) = l;
            L(i+1,:) = -l;
            i = i+2;
            L_face(f, :) = normal*K;
            
            % Collocation accross face direction ij
            cpt = G.cells.centroids(c1, :);  %centerpoint of cell 1
            
            %isJump = jumpFlag(selfFaces);
            %fa = selfFaces(isJump);
            %cells = N(selfFaces(~isJump), :);
            %cells(cells == c) = 0;
            %cells = sum(cells, 2);
            
            [T, ~, coeff, ~, successCell] = collocateSetOPT(G,cpt,selfFaces1,hap(selfFaces1,:), l);
            
            assert(successCell)
            
            A_ij(iF,1:nF_c1) = coeff.*var.O_j(selfFaces1);
            L1_ij(iF,1:nF_c1) = A_ij(iF,1:nF_c1)'.*var.O_i(selfFaces1);
            L2_ij(iF,1:nF_c1) = A_ij(iF,1:nF_c1)'.*var.O_j(selfFaces1);
            
            
            for j = 1:nF_c1
                T_all{j, 1}(f, :) = T(j, :);
            end
            
            coeff_ij(iF,1:nF_c1) = coeff;
            n = N(selfFaces1,:); n(n==c1)=[]; %neighbors to cell, including "ouside" cells
            act_cell_ij(iF,find(n~=0)) = true; %active cell neighbors
            n(n==0)=1; %set outside cell to 1 to avoid nonlogical indexing
            cn_1(iF,1:nF_c1) = n;  
           
            % Collocation accross face direction ji
           
            cpt = G.cells.centroids(c2,:);
           
            [T, ~, coeff, ~, successFace] = collocateSetOPT(G, cpt, selfFaces2, hap(selfFaces2,:), -l);
            
            assert(successFace)
            
            A_ji(iF,1:nF_c2) = coeff.*var.O_j(selfFaces2);
            L1_ji(iF,1:nF_c2) = A_ji(iF,1:nF_c2)'.*var.O_i(selfFaces2);
            L2_ji(iF,1:nF_c2) = A_ji(iF,1:nF_c2)'.*var.O_j(selfFaces2);
            
            for j = 1:nF_c2
                T_all{j, 2}(f, :) = T(j, :);
            end
            coeff_ji(iF,1:nF_c2) = coeff;
            n = N(selfFaces2,:); n(n==c2) = [];
            act_cell_ji(iF,find(n~=0)) = true;
            n(n==0)=1; %set outside cell to 1 to avoid nonlogical indexing
            cn_2(iF, 1:nF_c2) = n; 
        
        end
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
    coSet.C = {coeff_ij, coeff_ji};
   % coSet.types = {type_ij, type_ji};
    coSet.faceCellNeighbors = {cn_1, cn_2, cn_bf};
    coSet.l = L;
    coSet.active = struct('act_cell_ij',act_cell_ij,'act_cell_ji',act_cell_ji,...
        'act_cell_nghbr_i',act_cell_nghbr_i,'act_cell_nghbr_j',act_cell_nghbr_j,'act_cell_bf',act_cell_bf);
    
    
%     exclude_cell = struct('ix', faceNo, 'type', 2);
%     exclude_face = struct('ix', cellNo, 'type', 1);

    % Construct mapping operators
%     ops = cell(dim, 2);
%     Pc = speye(G.cells.num);
%     Pf = getFaceFromCellInterpolator(G, rock);
%     Pn = getNodeFromCellInterpolator(G);
    
%     for i = 1:numFace
%         for j = 1:2
%             if j == 1
%                 exclude = G.faces.neighbors(faceNo, 2);
%             else
%                 exclude = G.faces.neighbors(faceNo, 1);
%             end
%             gi = coSet.globalIndices{j}(:, i);
%             ti = coSet.types{j}(:, i);
%             ops{i, j} = getInterpolationOperatorOPT(Pc, Pf, Pn, gi, ti, exclude);
%         end
%     end
%     coSet.pressureOperators = ops;
    
    
   % coSet = storeFaceSet(G, rock, coSet, faceSign, L_face);
    
%     
%     intx = all(G.faces.neighbors > 0, 2);
%     ops = cell(dim, 2);
%     for i = 1:dim
%         for j = 1:2
%             if j == 1
%                 exclude = G.faces.neighbors(intx, 2);
%             else
%                 exclude = G.faces.neighbors(intx, 1);
%             end
%             gi = coSet.faceSet.globalIndices{j}(:, i);
%             ti =  coSet.faceSet.types{j}(:, i);
%             ops{i, j} = getInterpolationOperator(Pc, Pf, Pn, gi, ti, exclude);
%         end
%     end
%     coSet.faceSet.pressureOperators = ops;
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

function [T, types, coefficients, globIx, ok] = collocateSet(G, centerpt, cells, faces, nodes, l)
    [pts, typ, glob] = getCandidatePoints(G, cells, faces, nodes);
    [T, ix, coefficients, ok] = getSetForHF(centerpt, l, pts);
    globIx = glob(ix);
    types = typ(ix);
end

function [T, types, coefficients, globIx, ok] = collocateSetOPT(G,cpt,faces,hap, l)
   % [pts, typ, glob] = getCandidatePointsOPT(G,faces);
    warning('off','MATLAB:singularMatrix');
    [T, ix, coefficients, ok] = OPTgetSetForHF(cpt, l, hap);
    globIx = faces(ix); 
    %types = typ(ix);
    types = ones(4,1)*2;
end
 
% function [alpha,gamma] = decomp(faceNormal,noFaces,centroid,hap)
% 
% alpha = zeros(noFaces+1,1);
% sigma = zeros(noFaces,1);
% t = bsxfun(@minus,hap,centroid)'; %vec with x_k-x_i
% for i=1:length(t)
%     sigma(i) = norm(t(:,i))/norm(faceNormal);
%     t(:,i) = t(:,i)/norm(t(:,i));
% end
% t(:,noFaces+1)=0;
% fn = faceNormal/norm(faceNormal);
% 
% k = noFaces/min(sigma)+2; 
% c = ones(noFaces+1,1);
% c(end) = k;
% funk = @(x) c'*x;
% %funk = @(x) sum(x);
% option = optimset('Display','off');
% temp = ones(noFaces+1,1); temp(end)=0;
% lb = @(g) -g*temp; 
% alpha = fmincon(funk,alpha,[],[],t, fn',lb(alpha(end)),Inf*ones(noFaces+1,1),[],option);
% gamma = alpha(end);
% alpha(end)=[];
% alpha = alpha./sigma;
% 
% 
% %d = t*alpha;
% 
% end

function [T, ixSet, coeff, done] = OPTgetSetForHF(center,l,hap)
done = false;
ixSet = (1:1:length(hap));
%dim = size(pts,2);
N = size(hap,1);  %numper of faces to cell
coeff = zeros(N+1,1);

t = bsxfun(@minus,hap,center);
t_v = sqrt(sum(t.^2,2));   %norm of each row in t
t_u = bsxfun(@rdivide, t, t_v)';  %normalized t

l_v = norm(l,2);  %norm of 'facenormal'
l_u = l./l_v;     %normalized 'facenormal'

s = bsxfun(@rdivide,t_v,l_v);  %sigma
k = N/min(s)+1;              %kappa
c = ones(N+1,1);
c(end) = k;

%     if dim == 3
%     	nComb = (N-1)*(N-2)*N;
%     else
%         nComb = (N-1)*N;
%     end

% xl = center + l_u;  %???
% xi = bsxfun(@plus, center, t_u');    %???
% 
% dist = bsxfun(@minus, xl, xi);
% %     dist = cross(xi, repmat(xl, size(xi, 1), 1));
% dist = sum(dist.^2, 2);
% dist = sqrt(sum(dist.^2, 2));

%[~, ixSet] = sort(dist);



funk = @(x) c'*x;
temp = ones(N+1,1); temp(end)=0;
lb = @(g) -g*temp;
t_u(:,N+1)=0;

%if dim == 2  same for 3D ? 
option = optimset('Display','off');
coeff = fmincon(funk,coeff,[],[], t_u, l_u',lb(coeff(end)),Inf*ones(N+1,1),[],option);
gamma = coeff(end);
coeff(end) = [];
coeff = coeff./s;
if all(coeff >= -gamma)
    done = true;
end
%funk = @obj;
%lineq = struct('A', t,'b', faceNormal');
%[~,alpha,~] = unitBoxBFGS(alpha, funk, 'lineq', lineq, 'plotEvolution', false);
%else
%end

T = t(ixSet',:);

% for j=1:length(ixSet)
%     if coeff(j) == 0
%         ixSet(ixSet == j) = [];  %unnecessary? keep order
%     end
% end

end

function [f, df] = obj(x,c)
    x = initVariablesADI(x);
    f = c'*x;
    df = full(f.jac{1})';
    f = double(f);
end



function [T, ixSet, coeff, done] = getSetForHF(center, l, pts)
    dim = size(pts, 2);
%     assert(dim == 3);
    N = size(pts, 1);
    
    t = bsxfun(@minus, pts, center);
    t_v = sqrt(sum(t.^2, 2));
    t_u = bsxfun(@rdivide, t, t_v);
    
    l_v = norm(l, 2);
    l_u = l./l_v;
    
    xl = center + l_u;
    xi = bsxfun(@plus, center, t_u);
    
    
    dist = bsxfun(@minus, xl, xi);
%     dist = cross(xi, repmat(xl, size(xi, 1), 1));
    dist = sum(dist.^2, 2);
    dist = sqrt(sum(dist.^2, 2));
    
    [v, ix] = sort(dist);
    
    
    done = false;
    ptr = 1;
    
    if dim == 3
    	nComb = (N-1)*(N-2)*N;
    else
        nComb = (N-1)*N;
    end
    [coefficents, allsets] = deal(nan(nComb, dim));
    if dim == 3
        for i = 1:(N-2)
            if done
                break;
            end
            for j = (i+1):(N-1)
                if done
                    break;
                end
                for k = (j+1):N
                    I = ix(i);
                    J = ix(j);
                    K = ix(k);

                    t1 = t_u(I, :);
                    t2 = t_u(J, :);
                    t3 = t_u(K, :);

                    [C, D] = getCoefficients3D(t1, t2, t3, l_u);
                    coefficents(ptr, :) = C;

                    allsets(ptr, :) = [I, J, K];

                    if all(C >= 0) %&& abs(D) > 1e-8
                        if all(C <= 1)
                            done = true;
                            break
                        end
                    else
                        coefficents(ptr, :) = inf;
                    end
                    ptr = ptr + 1;
                end
            end
        end
    else
        for i = 1:(N-1)
            if done
                break;
            end
            for j = (i+1):N
                if done
                    break;
                end
                I = ix(i);
                J = ix(j);
                t1 = t_u(I, :);
                t2 = t_u(J, :);
                [C, D] = getCoefficients2D(t1, t2, l_u);
                coefficents(ptr, :) = C;
                allsets(ptr, :) = [I, J];
                if all(C >= 0)% && abs(D) > 1e-8
                    if all(C <= 1)
                        done = true;
                        break
                    end
                else
                    coefficents(ptr, :) = inf;
                end
                ptr = ptr + 1;
            end
        end
    end


    
    if done
        ixSet = allsets(ptr, :);
        coeff = coefficents(ptr, :);
    else
%         assert(ptr == nComb) 
        % Use fallback 
        mx = max(coefficents, [], 2);
        [v, sb] = min(mx);
        coeff = coefficents(sb, :);
        ixSet = allsets(sb, :);
        done = all(isfinite(coeff));
    end

    T = t(ixSet', :);
    for i = 1:numel(coeff)
        coeff(i) = coeff(i)./norm(T(i, :), 2);
    end
%     coeff = bsxfun(@rdivide, coeff, sqrt(sum(T.^2, 2)));;
%     T = bsxfun(@rdivide, T, sqrt(sum(T.^2, 2)));
    assert(all(coeff >= 0));
end

function D = computeDCoefficient3D(a, b, c)
    D = (dot(cross(a, b), c))/(norm(a, 2)*norm(b, 2)*norm(c, 2));
end

function D = computeDCoefficient2D(a, b)
    D = norm(a - b);
%     D = cross(a, b)
%     D = norm(a - b, 2)/(norm(a,2)*norm(b, 2));
end

function [C, D] = getCoefficients2D(t1, t2, l)
    D = nan;
    C = [0, 0];
    
    t1 = [t1, 0];%./norm(t1, 2);
    t2 = [t2, 0];%./norm(t2, 2);
    l = [l, 0];%./norm(l, 2);
    n = [0, 0, 1];
    C(1) = dot(n, cross(t2, l))/dot(n, cross(t2, t1));
    C(2) = dot(n, cross(t1, l))/dot(n, cross(t1, t2));
%     C = round(C, 8);
%     
%     l_new = t1*C(1) + t2*C(2);
%     C = C./norm(l_new, 2);
    l_new = t1*C(1) + t2*C(2);
%     assert(norm(l - l_new)/norm(l) < 1e-12);
%     C
%     D1 = computeDCoefficient2D(l, t2);
%     D2 = computeDCoefficient2D(t1, l);
% 
%     D = computeDCoefficient2D(t1, t2);
%     
%     C = [D1, D2]./D;
end


function [C, D] = getCoefficients3D(t1, t2, t3, l)
    D1 = computeDCoefficient3D(l, t2, t3);
    D2 = computeDCoefficient3D(t1, l, t3);
    D3 = computeDCoefficient3D(t1, t2, l);

    D = computeDCoefficient3D(t1, t2, t3);
    
    C = [D1, D2, D3]./D;
end


function [pts, types, ix] = getCandidatePoints(G, cells, faces, nodes)
    ix = [cells; faces; nodes];
    ix = ix(:);
    pts = [G.cells.centroids(cells, :); ...
           G.faces.centroids(faces, :); ...
           G.nodes.coords(nodes, :)];
    types = rldecode((1:3)', [numel(cells); numel(faces); numel(nodes)]);
end

function [pts, types, ix] = getCandidatePointsOPT(G,faces)
    ix = faces; %=globIx  - forenkle! 
    ix = ix(:);
    pts = G.faces.centroids(faces,:);
    types = ones(length(ix),1)*2;
end

function coSet = storeFaceSet(G, rock, coSet, faceSign, L_face)
    dim = G.griddim;
    [C_l, C_r, types_l, types_r, glob_l, glob_r, active_l, active_r] = ...
                                            deal(zeros(G.faces.num, dim));
    left = faceSign == 1;
    right = ~left;
    
    lf = G.cells.faces(left);
    rf = G.cells.faces(right);
    
    C_l(lf, :) = coSet.C{1}(left, :);
    C_r(rf, :) = coSet.C{1}(right, :);

    types_l(lf, :) = coSet.types{1}(left, :);
    types_r(rf, :) = coSet.types{1}(right, :);
    
    glob_l(lf, :) = coSet.globalIndices{1}(left, :);
    glob_r(rf, :) = coSet.globalIndices{1}(right, :);
    
    active_l(lf, :) = coSet.active{1}(left, :);
    active_r(rf, :) = coSet.active{1}(right, :);
    
    
    
    faceSet = struct();
    faceSet.C = {C_l, C_r};
    faceSet.types = {types_l, types_r};
    faceSet.globalIndices = {glob_l, glob_r};
    faceSet.active = {active_l, active_r};
    
    % Remove external
    intx = all(G.faces.neighbors > 0, 2);
    for i = 1:2
        faceSet.C{i} = faceSet.C{i}(intx, :);
        faceSet.types{i} = faceSet.types{i}(intx, :);
        faceSet.globalIndices{i} = faceSet.globalIndices{i}(intx, :);
        faceSet.active{i} = faceSet.active{i}(intx, :);
    end
    faceSet.l = L_face(intx, :);

    exclude_left = struct('ix', G.faces.neighbors(intx, 2), 'type', 1);
    exclude_right = struct('ix', G.faces.neighbors(intx, 2), 'type', 1);

    coSet.faceSet = faceSet;
end

function [hap,beta_ij,beta_ji,O_i,O_j] = findVariablesNTPFA(G,rock)

dim = G.griddim;
noFaces = G.faces.num;
bFace = boundaryFaces(G);
intFace = (1:1:noFaces)';
intFace(bFace) = [];
n_if = size(intFace,1); %number of internal faces

hap = zeros(noFaces,2);
hap(bFace,:) = G.faces.centroids(bFace,:); %hap for bf is equal to centroid

N = G.faces.neighbors;
nodePos = G.faces.nodePos;
nodes = G.faces.nodes;

beta_ij = zeros(noFaces,1); %beta for internal faces (left zero for boundary)
beta_ji = zeros(noFaces,1);
O_i = zeros(noFaces,1); %omega_ij (left zero for boundary)
O_j = zeros(noFaces,1); %omega_ji

dist = @(p1,p2,x) abs((p2(2)-p2(2))*x(1)-(p2(1)-p1(1))*x(2)+p2(1)*p1(2)-p2(2)*p1(1))/...
        sqrt((p2(2)-p1(2))^2+(p2(1)-p1(1))^2); 

for i = 1:n_if  %over all internal faces
    f = intFace(i);
    c1 = N(f,1); %cell neighbor 1
    c2 = N(f,2); %cell neighbor 2
    
    node = nodes(nodePos(f) : nodePos(f+1)-1, :);
    
    fn1 = G.nodes.coords(node(1),:);  %coordinates to node 1
    fn2 = G.nodes.coords(node(2),:);  %coordinates to node 2
    
    % for boundary faces, permeability is mirrored outside boundary
%     if c1 == 0 
%         K2 = getK(rock,c2,dim);
%         K1 = K2;
%         dist_j = abs((fn2(2)-fn1(2))*centroid2(1)-(fn2(1)-fn1(1))*centroid2(2)+fn2(1)*fn1(2)-fn2(2)*fn1(1))/...
%             sqrt((fn2(2)-fn1(2))^2+(fn2(1)-fn1(1))^2);
%         beta_ji(f) = K2*(G.faces.normals(f,:)*G.faces.normals(f,:)')/dist_j;
%     elseif c2 == 0
%         K1 = getK(rock,c1,dim);
%         K2 = K1;
%         dist_i = abs((fn2(2)-fn1(2))*centroid1(1)-(fn2(1)-fn1(1))*centroid1(2)+fn2(1)*fn1(2)-fn2(2)*fn1(1))/...
%             sqrt((fn2(2)-fn1(2))^2+(fn2(1)-fn1(1))^2);
%         beta_ij(f) = K1*(G.faces.normals(f,:)*G.faces.normals(f,:)')/dist_i;
%     else
%         K1 = getK(rock,c1,dim);
%         K2 = getK(rock,c2,dim);
%         centroid1 = G.cells.centroids(c1,:); %centroid of neighboring cell 1
%         centroid2 = G.cells.centroids(c2,:); %centroid of neighboring cell 2
%         dist_i = abs((fn2(2)-fn1(2))*centroid1(1)-(fn2(1)-fn1(1))*centroid1(2)+fn2(1)*fn1(2)-fn2(2)*fn1(1))/...
%             sqrt((fn2(2)-fn1(2))^2+(fn2(1)-fn1(1))^2);    % shortest distance from cell centroid 1 to face
%         dist_j = abs((fn2(2)-fn1(2))*centroid2(1)-(fn2(1)-fn1(1))*centroid2(2)+fn2(1)*fn1(2)-fn2(2)*fn1(1))/...
%             sqrt((fn2(2)-fn1(2))^2+(fn2(1)-fn1(1))^2);    % shortest distance from cell centroid 2 to face
%         beta_ij(f) = K1*(G.faces.normals(f,:)*G.faces.normals(f,:)')/dist_i;
%         %faceAreas(i)^2*
%         beta_ji(f) = K2*(G.faces.normals(f,:)*G.faces.normals(f,:)')/dist_j;
%         %faceAreas(i)^2*
%         O_i(f) = beta_ij(f)/(beta_ij(f)+beta_ji(f));  %omega
%         O_j(f) = beta_ji(f)/(beta_ij(f)+beta_ji(f));
%         hap(f,:) = (beta_ij(f)+beta_ji(f))\...
%             (beta_ij(f)*centroid1'+beta_ji(f)*centroid2'+(K1-K2)*G.faces.normals(f,:)');
%         
%     end
    
    % Find permeability
    K1 = getK(rock,c1,dim);
    K2 = getK(rock,c2,dim);
 
    % Centroids of neighboring cells:
    centroid1 = G.cells.centroids(c1,:); 
    centroid2 = G.cells.centroids(c2,:); 
    
    % Shortest distanse from cell centroids to face:
    dist_i = dist(fn1,fn2,centroid1);  
    dist_j = dist(fn1,fn2,centroid2);    
    
    % Calculate beta
    beta_ij(f) = K1*(G.faces.normals(f,:)*G.faces.normals(f,:)')/dist_i;
    %faceAreas(i)^2*
    if beta_ij(f) == inf
        beta_ij(f) = 0;
    end
    beta_ji(f) = K2*(G.faces.normals(f,:)*G.faces.normals(f,:)')/dist_j;
    %faceAreas(i)^2*
    if beta_ji(f) == inf
        beta_ji(f) = 0;
    end
    
    % Calculate omega:
    O_i(f) = beta_ij(f)/(beta_ij(f)+beta_ji(f));  
    O_j(f) = beta_ji(f)/(beta_ij(f)+beta_ji(f));
    
    % Calculate hap:   %HAP OUTSIDE GRID FOR SEAMOUNT!?!
    hap(f,:) = (beta_ij(f)+beta_ji(f))\...
        (beta_ij(f)*centroid1'+beta_ji(f)*centroid2'+(K1-K2)*G.faces.normals(f,:)');
    if any(isnan(hap(f,:)))
        hap(f,:) = G.faces.centroids(f,:);
    end
        
end


end