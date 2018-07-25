%% Example for NTPFAopt. 
%  We construct a twisted cartesian grid and consider a linear pressure 
%  drop from left to right. 
%  The NTPFAopt method is applied on the grid, and the solution is
%  shown. Also plotted is the pressure distribution along the x-axis.

%{
 Written by Pia-Kristina Heigrestad
%}

mrstModule add NTPFAopt ad-core solvers;

% Grid: 

n = 10;
G = computeGeometry(cartGrid([n,n], [1,1]));
G = computeGeometry(twister(G));

% Make rock and fluid:

rock = makeRock(G, 100*milli*darcy,0.002);
fluid = initSimpleADIFluid();

% Add boundary conditions:

[bc,src,W]=deal([]);

bc = pside(bc, G, 'west', 200*barsa, 'sat', [1,0]);
bc = pside(bc, G, 'east', 0, 'sat', [1,0]);

% Assemble model and solve problem: 

state0 = initResSol(G, 0, [0 1]);
model = PressureOilWaterModelNTPFAopt(G,rock,fluid);
state = incompSinglePhaseNTPFA(model, state0,'bc', bc, 'src',src);

% Plot solution: 

figure(1)
clf
plotCellData(G,state.pressure)
title('NTPFAopt')
axis equal tight
colorbar('Location','Southoutside')

figure(2)
clf
x = G.cells.centroids(:, 1);
plot(x,state.pressure,'.');
title('Pressure distribution along x-axis')
xlabel('cell position along x-axis')
ylabel('pressure')


