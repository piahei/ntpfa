dims = [5, 5];
pdims = [1, 1];
G = cartGrid(dims,pdims);  %create simple cartesian grid
G = computeGeometry(G);  %computes cell centroids and volumes, and face areas, centroids and normals.
rock = makeRock(G, 100*milli*darcy, 0.002);  %add rock properties
fluid = initSimpleADIFluid();  %initiate fluid
[bc, src, W] = deal([]);  %initiate boundary conditions, source, and wells
bc = pside(bc, G, 'Left', 1, 'sat', [1]);
bc = pside(bc, G, 'Right', 0, 'sat', [1]);
src = addSource(src, [1,G.cells.num], [1,1], 'sat', [2]);

state0 = initResSol(G, 0, [0 1]); %initial value
model = PressureOilWaterModelNTPFAopt(G,rock,fluid);
state = incompSinglePhaseNTPFA(model, state0,'bc', bc, 'src',src,'wells',W);

plotCellData(G,state.pressure)
colorbar('Location','eastoutside')