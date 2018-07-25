function [T, flux] = computeNonLinearTransForOpt(G, coSet, cellPressure,N)
% compute nonlinear transmissibilities for NTPFA using optimization
% in preprocessing step. Function is called in PressureOilWaterModelNTPFAopt
% SYNOPSIS:
%   [T,flux] = computeNonLinearTransForOpt(G, coSet, cellPressure, N)
%
% PARAMETERS:
%   G            - Grid structure as described by grid_structure.
%   coSet        - collocation set as produced by getCollactionSetOPT
%   cellPressure - ADI of lengt equal to number of cells with
%                  pressure values
%   N            - number-of-internal-faces x 2 vector with neighbors to
%                  each internal face
%
% RETURNS:
%   T    - structure with transmissibilities both ways across internal faces
%   flux - flux given as NaN
%
% SEE ALSO:
%   'getCollactionSetOPT', 'PressureOilWaterModelNTPFAopt'.

%{
Written by Pia-Kristina Heigrestad, 2018 
%}

if min(double(cellPressure)) < 0
    warning('Negative pressure in cells. Will fall back to linear TPFA.');
end
flux = nan;

T_if = computeTransIntFaces(coSet,cellPressure,N);
%T_bf = computeTransBoundaryFaces(G,coSet,cellPressure);

T = {T_if{1},T_if{2}};

end

function T_face = computeTransIntFaces(coSet,p,N)

c_i = N(:,1);
c_j = N(:,2);

cn_i = coSet.cellNeighbors{1};
cn_j = coSet.cellNeighbors{2};
act_ij = coSet.active.act_ij;
act_ji = coSet.active.act_ji;
A_ij = coSet.variables.A_ij;
A_ji = coSet.variables.A_ji;

pressure_i = p.val(cn_i).*act_ij; %set non-active cells to false
pressure_j = p.val(cn_j).*act_ji; %set non-active cells to false

L_ij = double2ADI(sum(A_ij.*pressure_i,2),p);
L_ji = double2ADI(sum(A_ji.*pressure_j,2),p);

% Determine weights:

[mu_ij, mu_ji] = findWeights(L_ij.val,L_ji.val);

%debugging
if any (mu_ij+mu_ji) ~= 1
    warning('fault in weights')
end

% Determine eta:

e_ij = mu_ij.*sum(A_ij,2) + mu_ji.*sum(A_ji.*~act_ji,2);
e_ji = mu_ji.*sum(A_ji,2) + mu_ij.*sum(A_ij.*~act_ij,2);

% Determine residual:
r = mu_ji.*L_ji - mu_ij.*L_ij;

a_r = abs(r);

% Assemble transmissibility:

T_1 = e_ij+(a_r+r)./((p(c_i)+eps).*2);
T_2 = e_ji+(a_r-r)./((p(c_j)+eps).*2);

T_face = {T_1,T_2};

end

function [w_ij,w_ji] = findWeights(L_ij,L_ji)

[w_ij,w_ji] = deal(zeros(numel(L_ij),1));

nonzeroRow = any([L_ij,L_ji],2);
zeroRow = ~nonzeroRow;
nonzeroRow = find(nonzeroRow);
zeroRow = find(zeroRow);
%absolute value for weights
a_ij = abs(L_ij);
a_ji = abs(L_ji);

% if L_ij=L_ji=0 then mu_ij=mu_ji=0.5:

w_ij(nonzeroRow) = a_ji(nonzeroRow)./(a_ij(nonzeroRow)+a_ji(nonzeroRow));
w_ij(zeroRow) = 0.5;
w_ji(nonzeroRow) = a_ij(nonzeroRow)./(a_ij(nonzeroRow)+a_ji(nonzeroRow));
w_ji(zeroRow) = 0.5;

end

% function T_bf = computeTransBoundaryFaces(G,coSet,p)
% 
% L1 = coSet.variables.boundaryVariables.L1_bf;
% L2 = coSet.variables.boundaryVariables.L2_bf;
% A = coSet.variables.boundaryVariables.A_bf;
% act_cells = coSet.active.act_cell_bf;
% 
% isBF = false(G.faces.num, 1);
% isBF(boundaryFaces(G)) = true;
% nghbrs = coSet.cellNeighbors{3};
% pressure = p.val(nghbrs);
% pressure = bsxfun(@times,pressure,act_cells);
% c = sum(G.faces.neighbors(isBF,:),2);
% 
% L = sum(bsxfun(@times,L1,p.val(c)),2)+sum(bsxfun(@times,L2,pressure),2);
% 
% e = sum(A,2);
% 
% T_bf = e+(abs(L)+L)./(2*(p(c)+eps));
% 
% end
