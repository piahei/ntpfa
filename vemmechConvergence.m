%% Unit Square test on VEMMECH grid

d = [10 20 30 40];

fluid = initSimpleADIFluid();
fluid2 = initSimpleFluid('mu', [1, 1], 'rho', [1, 1], 'n', [1, 1]);

[error_opt,error_lin,error_mpfa] = deal(zeros(length(d),1));

p_ex = @(cntr) sin(cntr(1))+2*cntr(2);

for i=1:length(d)
    opt = struct('L'         , [d(i),d(i)], ...
        'cartDims'  , [10 10], ...
        'grid_type' , 'square', ...
        'disturb'   , 0.1, ... %parameter for disturbing grid
        'E'         , 4e7, ...  %youngs modolo
        'nu'        , 0.44);% poiso ratio
    
    G = squareGrid(opt.cartDims,opt.L,'grid_type','mixed3','disturb',opt.disturb);
    G = computeGeometry(G);
    
    rock = makeRock(G, 100*milli*darcy, .02);
    %K = [1, 0.5; 0.5, 1];
    
    %rock.perm = repmat([K(1,1), K(1,2), K(2,2)], G.cells.num, 1);
    
    p_exact = zeros(G.cells.num,1);
    for j=1:G.cells.num
        p_exact(j) = p_ex(G.cells.centroids(j,:));
    end
    
    [Lx, Ly] = deal(opt.L(1), opt.L(2));
    assert(G.griddim == 2);
    x = [0, Lx];
    
    [bc,src,W]=deal([]);
    
    bndryFaces = boundaryFaces(G);
    for k = 1:length(bndryFaces)
        pres = p_ex(G.faces.centroids(k,:));
        bc = addBC(bc,bndryFaces(k),'pressure',pres,'sat',[1]);
    end
    
    state0 = initResSol(G, 0, [0 1]);
    
    model = PressureOilWaterModelNTPFAopt(G,rock,fluid);%,'Algorithm','interior-point'); %can use sqp, faster but worse result in some cases
    model2 = PressureOilWaterModelNTPFA(G,rock,fluid);
    modelMPFA = computeMultiPointTrans(G, rock);
    
    state = incompSinglePhaseNTPFA(model, state0,'bc', bc, 'src',src);
    state2 = incompSinglePhaseNTPFA(model2,state0,'bc',bc,'src',src);
    mpfa = incompMPFA(state0, G, modelMPFA, fluid2, 'bc', bc,'src', src);
    
    error_opt(i) = sqrt((sum(G.cells.volumes.*(p_exact-state.pressure).^2))/(sum(G.cells.volumes.*p_exact.^2)));
    error_lin(i) = sqrt((sum(G.cells.volumes.*(p_exact-state2.pressure).^2))/(sum(G.cells.volumes.*p_exact.^2)));
    error_mpfa(i) = sqrt((sum(G.cells.volumes.*(p_exact-mpfa.pressure).^2))/(sum(G.cells.volumes.*p_exact.^2)));
    
    if i==1
        figure(1)
        clf
        plotCellData(G,p_exact)
        title('exact pressure')
        
        figure(2)
        clf
        subplot(1,3,1)
        plotCellData(G,state.pressure)
        axis equal tight
        title('NTPFAopt')
        colorbar('Location','southoutside')
        
        subplot(1,3,2)
        plotCellData(G,state2.pressure)
        axis equal tight
        title('NTPFA')
        colorbar('Location','southoutside')
        
        subplot(1,3,3)
        plotCellData(G,mpfa.pressure)
        axis equal tight
        title('MPFA')
        colorbar('Location','southoutside')
    end
end

nn = d.*d;
figure(3)
clf
loglog(nn,error_opt)
hold on
loglog(nn,error_lin)
loglog(nn,error_mpfa)
axis tight
legend('Error NTPFAopt','Error NTPFA','Error MPFA')
hold off

