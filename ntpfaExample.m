mrstModule add book pia ad-core ad-blackoil ad-props solvers...
               blackoil-sequential mimetic incomp mpfa streamlines mrst-gui...
               vemmech;
%% Standard cart grid ------------------------------------------------------

% n = 10;
% G = computeGeometry(cartGrid([n,n], [1,1]));
% G = twister(G);
% rock = makeRock(G, 100*milli*darcy, 0.002);
% %%Anisotropic perm:-----------------
% K = [1, 0; 0, 100]*100*milli*darcy;
% R = @(t) [cos(t), -sin(t); sin(t), cos(t)];
% t = pi/8;
% K = R(t)*K*R(t)';
% rock.perm = repmat([K(1,1), K(1,2), K(2,2)], G.cells.num, 1);
% %-----------------------------------
% [src,W,bc] = deal([]);
% bc = pside(bc, G, 'west', 150*barsa, 'sat', [1,0]);
% bc = pside(bc, G, 'east', 0*barsa, 'sat', [1,0]);

%% glued grid --------------------------------------------------------------
% 
% G1 = computeGeometry(cartGrid([1 2],[1,1]));
% G2 = computeGeometry(cartGrid([1,3],[1,1]));
% G = glue2DGrid(G1,translateGrid(G2,[1,0]));
% G = computeGeometry(G);
% rock = makeRock(G,100*milli*darcy, 0.002);
% [src,W,bc] = deal([]);
% 
% % bc = addBC(bc,[1,3],'pressure',0,'sat',[1]);
% % bc = addBC(bc,[13,16,19],'pressure',50*barsa,'sat',[1]);
% 
% x = [0, 2];
% for i = 1 : 2
%     faces = find(abs(G.faces.centroids(:, 1) - x(i)) < eps);
%     bc = addBC(bc, faces, 'pressure', 0,'sat',[1]);
% %     bc{i} = addBC([],faces,'pressure',0);
% %     bc{i} = rmfield(bc{i}, 'type');
% %     bc{i} = rmfield(bc{i}, 'sat');
% end
% y = [0, 1];
% for i = 1 : 2
%     faces = find(abs(G.faces.centroids(:, 2) - y(i)) < eps);
%     bc = addBC(bc,faces,'pressure',100*barsa,'sat',[1]);
% %     bc{i + 2} = addBC([], faces, 'pressure', 50*barsa);
% %     bc{i + 2} = rmfield(bc{i + 2}, 'type');
% %     bc{i + 2} = rmfield(bc{i + 2}, 'sat');
% end
% % --------------------------------------------------------------------------
%% Skew grid ---------------------------------------------------------------
% G = cartGrid([41,20],[2,1]);
% makeSkew = @(c) c(:,1) + .4*(1-(c(:,1)-1).^2).*(1-c(:,2));
% G.nodes.coords(:,1) = 2*makeSkew(G.nodes.coords);
% G = computeGeometry(G);
% rock = makeRock(G, 100*milli*darcy, 0.2);
% pv   = sum(poreVolume(G,rock));
% [bc, src, W] = deal([]);
% srcCells = findEnclosingCell(G,[2 .975; .5 .025; 3.5 .025]);
% src = addSource(src, srcCells, [1; -0.5; -0.5],'sat',[1]);

%% Vemmech grid -------------------------------------------------------------
opt = struct('L'         , [2 2], ...
             'cartDims'  , [4 4], ...
             'grid_type' , 'square', ...
             'disturb'   , 0.1, ... %parameter for disturbing grid
             'E'         , 4e7, ...  %youngs modolo
             'nu'        , 0.44);% poiso ratio

G = squareGrid(opt.cartDims,opt.L,'grid_type','mixed2','disturb',opt.disturb);
G = computeGeometry(G);

rock = makeRock(G, 100*milli*darcy,0.002);

[Lx, Ly] = deal(opt.L(1), opt.L(2));
assert(G.griddim == 2);
x = [0, Lx];

[bc,src,W]=deal([]);

% bc = cell(4,1);
for i = 1 : 2
    faces = find(abs(G.faces.centroids(:, 1) - x(i)) < eps);
    bc = addBC(bc, faces, 'pressure', 0,'sat',[1]);
%     bc{i} = addBC([],faces,'pressure',0);
%     bc{i} = rmfield(bc{i}, 'type');
%     bc{i} = rmfield(bc{i}, 'sat');
end
y = [0, Ly];
for i = 1 : 2
    faces = find(abs(G.faces.centroids(:, 2) - y(i)) < eps);
    bc = addBC(bc,faces,'pressure',50*barsa,'sat',[1]);
%     bc{i + 2} = addBC([], faces, 'pressure', 50*barsa);
%     bc{i + 2} = rmfield(bc{i + 2}, 'type');
%     bc{i + 2} = rmfield(bc{i + 2}, 'sat');
end

%---------------------------------------------------------------

%% Seamount -----------------------------------------------------------------
% load seamount
% %  
% T = triangleGrid([x(:) y(:)], delaunay(x,y));
% [Tmin,Tmax] = deal(min(T.nodes.coords), max(T.nodes.coords));
% T.nodes.coords = bsxfun(@times, ...
% bsxfun(@minus, T.nodes.coords, Tmin), 1000./(Tmax - Tmin));
% T = computeGeometry(T); 
% % 
% % % Triangular:
%  G = T;

% Cartesian:
% G = computeGeometry(cartGrid([25 25], [1000 1000]));
% inside = isPointInsideGrid(T, G.cells.centroids);
% G = removeCells(G, ~inside);

%Radial/PEBI: PROBLEMER HER! HAP = CENTROID->NaN
% P = [];
% for r = exp([-3.5:.2:0, 0, .1]),
% [x,y] = cylinder(r,25); P = [P [x(1,:); y(1,:)]];
% end
% P = unique([P'; 0 0],'rows');
% [Pmin,Pmax] = deal(min(P), max(P));
% P = bsxfun(@minus, bsxfun(@times, ...
% bsxfun(@minus, P, Pmin), 1200./(Pmax-Pmin)), [150 100]);
% inside = isPointInsideGrid(T, P);
% G = computeGeometry( pebi( triangleGrid(P(inside,:)) ));


% G = triangleGrid([x(:) y(:)], delaunay(x,y));
% [Gmin,Gmax] = deal(min(G.nodes.coords), max(G.nodes.coords));
% G.nodes.coords = bsxfun(@times, ...
% bsxfun(@minus, G.nodes.coords, Gmin), 1000./(Gmax - Gmin));
% G = computeGeometry(G);

% properties:
% rock = makeRock(G, 100*milli*darcy, 0.2);
% e=boundaryFaces(G);
% [bc,src,W]=deal([]);
% bc = addBC(bc, e, 'pressure', 50*barsa, 'sat',[1]);
% tmp = (G.cells.centroids - repmat([450, 500],G.cells.num,1)).^2;
% [~,ind] = min(sum(tmp,2));
% pv = sum(poreVolume(G,rock));
%  %src = addSource(src, ind, -.02*pv/year,'sat',[1]);
% src = addSource(src, ind, -.02*pv/year,'sat',[1]);


%  rock.poro = repmat(0.2, G.cells.num, 1);
%    rock.perm = repmat(100*milli*darcy, G.cells.num, 1);


%% ----------------------------------------------------------------


 %%SEAMOUNT--------------------------------------------------------
% load seamount
% G = triangleGrid([x(:) y(:)], delaunay(x,y));
% [Gmin,Gmax] = deal(min(G.nodes.coords), max(G.nodes.coords));
% G.nodes.coords = bsxfun(@times, ...
% bsxfun(@minus, G.nodes.coords, Gmin), 1000./(Gmax - Gmin));
% G = computeGeometry(G);
% rock = makeRock(G, 0.05, 0.2);
% e=boundaryFaces(G);
% [bc,src,W]=deal([]);
% bc = addBC(bc, e, 'pressure', 50*barsa, 'sat',[1]);
% tmp = (G.cells.centroids - repmat([450, 500],G.cells.num,1)).^2;
% [~,ind] = min(sum(tmp,2));
% %src = addSource(src, ind, -.02*pv/year,'sat',[1]);
% src = addSource(src, ind, -.02,'sat',[1]);
%-------------------------------------------------------------------


%% Isotropic perm

fluid = initSimpleADIFluid();
%  nkr = [1,1];
%  mu  = [0.1, 1]*centi*poise;
%  rho = [1000,800]*kilogram/meter^3;
%  fluid = initSimpleADIFluid('phases', 'WO', 'mu', mu, 'rho', rho, 'n', nkr);

%state0 = initResSol(G, 100*barsa, [0,1]);
state0 = initResSol(G, 0, [1 0]);

% [src,W,bc] = deal([]);
% bc = pside(bc, G, 'west', 100*barsa, 'sat', [1,0]);
% bc = pside(bc, G, 'east', 50*barsa, 'sat', [1,0]);

model = PressureOilWaterModelNTPFAopt(G,rock,fluid);  
model2 = PressureOilWaterModelNTPFA(G,rock,fluid);
%dT = 1;
%schedule = simpleSchedule(dT, 'W', W, 'bc', bc, 'src', src);
state = incompSinglePhaseNTPFA(model, state0,'bc', bc, 'src',src);
state2 = incompSinglePhaseNTPFA(model2,state0,'bc',bc,'src',src);


figure(1)
clf
subplot(1,2,1)

plotCellData(G,state.pressure)
% hold on
% plot([.5 2 3.5], [.025 .975 .025],'.','Color',[.9 .9 .9],'MarkerSize',16);
% hold off
title('NTPFA OPT')
axis equal
h=colorbar('Location','South');
subplot(1,2,2)
plotCellData(G,state2.pressure)
% hold on
% plot([.5 2 3.5], [.025 .975 .025],'.','Color',[.9 .9 .9],'MarkerSize',16);
% hold off
title('NTPFA lin')
axis equal
h=colorbar('Location','South');
%plotToolbar(G, state)


figure(2)
clf
x = G.cells.centroids(:,1);
plot(x,state.pressure,'.');
hold on
plot(x,state2.pressure,'.');
legend('ntpfaOPT','ntpfaLIN');

%% Anisotropic perm

% K = [1, 0; 0, 100]*100*milli*darcy;
% R = @(t) [cos(t), -sin(t); sin(t), cos(t)];
% t = pi/8;
% K = R(t)*K*R(t)';
% 
% rock.perm = repmat([K(1,1), K(1,2), K(2,2)], G.cells.num, 1);
% 
% model = PressureOilWaterModelNTPFA(G, rock, fluid);
% state = incompSinglePhaseNTPFA(model, state0, 'bc', bc);
% 
