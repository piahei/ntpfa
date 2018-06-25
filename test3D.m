%%3D example

% [nx,ny,nz] = deal(20, 20, 5);
% [Lx,Ly,Lz] = deal(1000, 1000, 50);

[nx,ny,nz] = deal(2, 2, 2);
[Lx,Ly,Lz] = deal(1, 1, 1);

G = cartGrid([nx ny nz], [Lx Ly Lz]);
G.nodes.coords(:,3) = G.nodes.coords(:,3) + 500;
G = computeGeometry(G);
bc=[];
% bc = fluxside(bc, G, 'EAST', 5e3*meter^3/day);
% bc = pside (bc, G, 'WEST', 50*barsa);

bc = pside(bc,G,'EAST', 0);
bc = pside(bc,G,'WEST', 50*barsa);

src = [];

rock = makeRock(G, 100*milli*darcy, .2);
fluid = initSimpleADIFluid();

state0 = initResSol(G, 0, [0 1]);
model = PressureOilWaterModelNTPFAopt(G,rock,fluid);
model2 = PressureOilWaterModelNTPFA(G,rock,fluid);

state = incompSinglePhaseNTPFA(model, state0,'bc', bc, 'src',src);
state2 = incompSinglePhaseNTPFA(model2,state0,'bc',bc,'src',src);

figure(1)
clf
plotCellData(G,state.pressure)
view(-125,20),camproj perspective
colorbar('Location','Southoutside');

figure(2)
clf
plotCellData(G,state2.pressure)
view(-125,20),camproj perspective
colorbar('Location','Southoutside');
