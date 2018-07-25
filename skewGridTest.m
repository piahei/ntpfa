%% Compare results for the TPFA and NTPFA schemes on a non K-orthogonal grid
% We consider a homogeneous, rectangular reservoir with a symmetric well
% pattern consisting of one injector and two producers. Because of the
% symmetry, the travel times from the injector to each producer should be
% equal. When using a skew grid that is not K-orhtogonal, the travel times
% of the TPFA method will not be equal and the flow pattern will differ 
% quite a lot from being symmetric. The NTPFAopt method produces symmetric,
% well defined results, whereas the NTPFA method yields inaccurate results.
mrstModule add incomp diagnostics streamlines pia ad-core

% Rectangular reservoir with a skew grid.
G = cartGrid([60+1,30],[2,1]);
makeSkew = @(c) c(:,1) + .4*(1-(c(:,1)-1).^2).*(1-c(:,2));
G.nodes.coords(:,1) = 2*makeSkew(G.nodes.coords);
G = computeGeometry(G);
disp(['  Grid: ' num2str(G.cartDims)]);

% Homogeneous reservoir properties
rock = makeRock(G, 100*milli*darcy, .2);
pv   = sum(poreVolume(G,rock));

% Symmetric well pattern
srcCells = findEnclosingCell(G,[2 .975; .5 .025; 3.5 .025]);
src = addSource([], srcCells, [pv; -.5*pv; -.5*pv],'sat',[1]);

% Single-phase fluid for TPFA
fluid = initSingleFluid('mu', 1*centi*poise,'rho', 1000*kilogram/meter^3);

% Solve flow problem using TPFA:
hT     = computeTrans(G, rock);
stateTPFA  = initState(G,[], 0);
stateTPFA  = incompTPFA(stateTPFA, G, hT, fluid, 'src', src);
tofTPFA    = computeTimeOfFlight(stateTPFA, G, rock, 'src', src);

% Solve flow problem NTPFAopt:
state0 = initResSol(G, 0, [1 0]);
fluidNTPFA = initSimpleADIFluid();
modelNTPFAopt = PressureOilWaterModelNTPFAopt(G,rock,fluidNTPFA);
stateNTPFAopt = incompSinglePhaseNTPFA(modelNTPFAopt, state0,'src',src);
tofNTPFAopt    = computeTimeOfFlight(stateNTPFAopt, G, rock, 'src', src);

% solve flow probem NTPFA:
modelNTPFA = PressureOilWaterModelNTPFA(G,rock,fluidNTPFA);
stateNTPFA = incompSinglePhaseNTPFA(modelNTPFA,state0,'src',src);
tofNTPFAlin = computeTimeOfFlight(stateNTPFA,G,rock,'src',src);

% Plot solution:
figure(1)
clf
subplot(3,1,1);
plotCellData(G,stateTPFA.pressure,'EdgeColor','k','EdgeAlpha',.05);
axis equal tight
hold on
plot([.5 2 3.5], [.025 .975 .025],'.','Color',[.9 .9 .9],'MarkerSize',16);
hold off

subplot(3,1,2);
plotCellData(G, stateNTPFAopt.pressure, 'EdgeColor', 'k', 'EdgeAlpha', .05);
axis equal tight
hold on
plot([.5 2 3.5], [.025 .975 .025],'.','Color',[.9 .9 .9],'MarkerSize',16);
hold off

subplot(3,1,3);
plotCellData(G,stateNTPFA.pressure,'EdgeColor','k','EdgeAlpha',.05);
axis equal tight
hold on
plot([.5 2 3.5], [.025 .975 .025],'.','Color',[.9 .9 .9],'MarkerSize',16);
hold off

figure(2)
subplot(3,1,1);
plotCellData(G, tofTPFA, tofTPFA<.2, 'EdgeColor','none'); caxis([0 .2]); box on
seed = floor(G.cells.num/5)+(1:G.cartDims(1))';
hf = streamline(pollock(G, stateTPFA, seed, 'substeps', 1) );
hb = streamline(pollock(G, stateTPFA, seed, 'substeps', 1, 'reverse' , true));
set ([ hf ; hb ],'Color','k');

subplot(3,1,2);
plotCellData(G, tofNTPFAopt, tofNTPFAopt<.2, 'EdgeColor','none'); caxis([0 .2]); box on
seed = floor(G.cells.num/5)+(1:G.cartDims(1))';
hf = streamline(pollock(G, stateNTPFAopt, seed, 'substeps', 1) );
hb = streamline(pollock(G, stateNTPFAopt, seed, 'substeps', 1, 'reverse' , true));
set ([ hf ; hb ], 'Color' , 'k' );

subplot(3,1,3);
plotCellData(G, tofNTPFAlin, tofNTPFAlin<.2, 'EdgeColor','none'); caxis([0 .2]); box on
seed = floor(G.cells.num/5)+(1:G.cartDims(1))';
hf = streamline(pollock(G, stateNTPFA, seed, 'substeps', 1) );
hb = streamline(pollock(G, stateNTPFA, seed, 'substeps', 1, 'reverse' , true));
set ([ hf ; hb ], 'Color' , 'k' );

%{
Copyright 2009-2017 SINTEF ICT, Applied Mathematics.

This file is part of The MATLAB Reservoir Simulation Toolbox (MRST).

MRST is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

MRST is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with MRST.  If not, see <http://www.gnu.org/licenses/>.
%}
