gravity off
mrstModule add ntpfa ad-core ad-blackoil ad-props ...
               blackoil-sequential mimetic incomp

dims = [15, 15];
pdims = [1, 1];

G = cartGrid(dims, pdims);
% G = twister(G, 0.1);
G = computeGeometry(G);
fluid = initSimpleADIFluid();

rock = makeRock(G, 1, 1);
model = SimplePressureModel(G, rock, fluid);
model.verbose = true;

[bc, src, W] = deal([]);
bc = pside(bc, G, 'xmin', 1, 'sat', [1]);
bc = pside(bc, G, 'xmax', 0, 'sat', [1]);

% Comment out for source
% bc = [];
% src = addSource(src, 1, 1, 'sat', [1]);
% src = addSource(src, G.cells.num, -1, 'sat', [1]);

%%
% NTPFA
state0 = initResSol(G, 0, [1, 0]);
dT = 1;

schedule = simpleSchedule(dT, 'W', W, 'bc', bc, 'src', src);
[~, states] = simulateScheduleAD(state0, model, schedule);
ntpfa = states{1};
% MIMETIC
S = computeMimeticIP(G, rock);
fluid2 = initSimpleFluid('mu', [1, 1], 'rho', [1, 1], 'n', [1, 1]);
mimetic = incompMimetic(state0, G, S, fluid2, 'bc', bc, 'src', src);

%%
figure(1), clf
subplot(1, 2, 1);
plotCellData(G, ntpfa.pressure)
title('TPFA')

subplot(1, 2, 2);
plotCellData(G, mimetic.pressure)
title('Mimetic')

