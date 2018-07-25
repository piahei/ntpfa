%% 3D convergence test

nx = [10,20,30,40];
ny = [10,15,20,25];
nz = [5,5,5,5];

[Lx,Ly,Lz] = deal(1, 1, 1);
x = [0, Lx];
y = [0, Ly];
z = [0, Lz];

error_opt = zeros(length(nx),1);
error_lin = zeros(length(nx),1);
error_mpfa = zeros(length(nx),1);

%%anisotropic p_ex:
p_ex = @(pnt) sin(2*pi*pnt(1))*sin(2*pi*(pnt(2)))*sin(2*pi*(pnt(3)))+1;
%%isotropic p_ex:
%p_ex = @(pnt) sin(pi*pnt(1))*sin(pi*(pnt(2)+0.5))*sin(pi*(pnt(3)+1/3))+1;
fluid = initSimpleADIFluid();
fluid2 = initSimpleFluid('mu', [1, 1], 'rho', [1, 1], 'n', [1, 1]);

for i = 1:length(nx)
    G = cartGrid([nx(i),ny(i),nz(i)], [Lx,Ly,Lz]);
    % G.nodes.coords(:,3) = G.nodes.coords(:,3) + 500;
    %G = twister(G);
    G = computeGeometry(G);
    
    p_exact = zeros(G.cells.num,1);
    for j=1:G.cells.num
        p_exact(j) = p_ex(G.cells.centroids(j,:));
    end
    
    [bc,src,w] = deal([]);
    
    rock = makeRock(G, 100*milli*darcy, .02);
    %K = [1, 0.5, 0; 0.5, 1, 0.5; 0, 0.5, 1]*milli*darcy;
    K = [1, 0, 0; 0, 1, 0; 0, 0, 10^3];%*milli*darcy;
    rock.perm = repmat([K(1,1), K(1,2), K(1,3),K(2,2),K(2,3),K(3,3)], G.cells.num, 1);
%     
    bndryFaces = boundaryFaces(G);
    for j = 1:length(bndryFaces);
        p_temp = p_ex(G.faces.centroids(bndryFaces(j),:));
        bc = addBC(bc, bndryFaces(j), 'pressure', p_temp,'sat',[1]);
    end

    state0 = initResSol(G, 0, [0 1]);
    model = PressureOilWaterModelNTPFAopt(G,rock,fluid);
    model2 = PressureOilWaterModelNTPFA(G,rock,fluid);
    modelMPFA = computeMultiPointTrans(G, rock);
    
    
    state = incompSinglePhaseNTPFA(model, state0,'bc', bc, 'src',src);
    state2 = incompSinglePhaseNTPFA(model2,state0,'bc',bc,'src',src);
    MPFA = incompMPFA(state0, G, modelMPFA, fluid2, 'bc', bc,'src', src);
    
    error_opt(i) = sqrt((sum(G.cells.volumes.*(p_exact-state.pressure).^2))/(sum(G.cells.volumes.*p_exact.^2)));
    error_lin(i) = sqrt((sum(G.cells.volumes.*(p_exact-state2.pressure).^2))/(sum(G.cells.volumes.*p_exact.^2)));
    error_mpfa(i) = sqrt((sum(G.cells.volumes.*(p_exact-MPFA.pressure).^2))/(sum(G.cells.volumes.*p_exact.^2)));
%     
    if i==2
        figure(1)
        clf
        subplot(1,4,1)
        plotCellData(G,p_exact)
        view(-124,20),camproj perspective
        axis equal tight
        title('Exact pressure')
        colorbar('Location','southoutside')
        
        subplot(1,4,2)
        plotCellData(G,state.pressure)
        view(-124,20),camproj perspective
        axis equal tight
        title('NTPFAopt')
        colorbar('Location','southoutside')
        
        subplot(1,4,3)
        plotCellData(G,state2.pressure)
        view(-124,20),camproj perspective
        axis equal tight
        title('NTPFAlin')
        colorbar('Location','southoutside')
        
        subplot(1,4,4)
        plotCellData(G,MPFA.pressure)
        view(-124,20),camproj perspective
        axis equal tight
        title('MPFA')   
        colorbar('Location','southoutside')
    end
 end

nn = nx.*ny.*nz;
figure(2)
clf
loglog(nn,error_opt)
hold on
loglog(nn,error_lin)
loglog(nn,error_mpfa)
%box off
axis tight
legend('error NTPFAopt','error NTPFAlin','error MPFA')




