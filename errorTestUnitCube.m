%% CONVERGENCE TEST ON UNIT SQUARE

n = [10,20,40,80];

error_opt = zeros(length(n),1);
error_lin = zeros(length(n),1);
%error_mpfa = zeros(length(n),1);

%p_ex = @(cntr) sin(pi*cntr(1))*sin(pi*cntr(2));
p_ex = @(cntr) sin(cntr(1))+2*cntr(2); 

%q = @(cntr) -sin(cntr(1))-0.5*sin(cntr(1));

fluid = initSimpleADIFluid();
fluid2 = initSimpleFluid('mu', [1, 1], 'rho', [1, 1], 'n', [1, 1]);

for i=1:length(n)
    G = computeGeometry(cartGrid([n(i),n(i)], [1,1]));
    G = computeGeometry(twister(G));
    
    p_exact = zeros(G.cells.num,1);
    for j=1:G.cells.num
        p_exact(j) = p_ex(G.cells.centroids(j,:));
    end
    
    rock = makeRock(G, 100*milli*darcy, 0.002);
    %%Anisotropic perm:-----------------
    K = [1, 0.5; 0.5, 1];
    % R = @(t) [cos(t), -sin(t); sin(t), cos(t)];
    % t = pi/8;
    % K = R(t)*K*R(t)';
    rock.perm = repmat([K(1,1), K(1,2), K(2,2)], G.cells.num, 1);
    
    %-----------------------------------
    [src,W,bc] = deal([]);
   
    bndryFaces=boundaryFaces(G);
    cntBndryFaces = G.faces.centroids(bndryFaces,:);
    for k = 1 : length(bndryFaces)
        p_ex_temp = p_ex(cntBndryFaces(k,:));
        bc = addBC(bc, bndryFaces(k), 'pressure', p_ex_temp,'sat',[1]);
    end
    
    state0 = initResSol(G, 0, [0 1]);
    model = PressureOilWaterModelNTPFAopt(G,rock,fluid);
    model2 = PressureOilWaterModelNTPFA(G,rock,fluid);
   % modelMPFA = computeMultiPointTrans(G, rock);
   % MPFA = incompMPFA(state0, G, modelMPFA, fluid2, 'bc', bc,'src', src);

    state = incompSinglePhaseNTPFA(model, state0,'bc', bc, 'src',src);
    state2 = incompSinglePhaseNTPFA(model2,state0,'bc',bc, 'src',src);
    
    error_lin(i) = sqrt((sum(G.cells.volumes.*(p_exact-state2.pressure).^2))/(sum(G.cells.volumes.*p_exact.^2)));
    error_opt(i) = sqrt((sum(G.cells.volumes.*(p_exact-state.pressure).^2))/(sum(G.cells.volumes.*p_exact.^2)));
   % error_mpfa(i) = sqrt((sum(G.cells.volumes.*(p_exact-MPFA.pressure).^2))/(sum(G.cells.volumes.*p_exact.^2)));
    if i==1
        figure(1)
        clf
        subplot(1,3,1)
        plotCellData(G,p_exact)
        title('exact')
        axis equal tight
        colorbar('location', 'southoutside');
        
        subplot(1,3,2)
        plotCellData(G,state.pressure)
        title('NTPFAopt')
        axis equal tight
        colorbar('Location','Southoutside');
        
        subplot(1,3,3)
        plotCellData(G,state2.pressure)
        title('NTPFAlin pressure')
        axis equal tight
        colorbar('Location','Southoutside');
        
%         subplot(1,4,4)
%         plotCellData(G,MPFA.pressure)
%         title('MPFA')
%         axis equal tight
%         colorbar('Location','Southoutside');
        
        figure(2)
        clf
        x = G.cells.centroids(:, 1);
        plot(x,state.pressure,'.');
        hold on
        plot(x,state2.pressure,'.');
        %plot(x,MPFA.pressure,'.');
        plot(x,p_exact,'.');
        legend('NTPFAopt','NTPFAlin','exact');
        hold off
    end
    
end

nn = n.*n;
figure(3)
clf
loglog(nn,error_lin)
hold on
loglog(nn,error_opt)
%loglog(nn,error_mpfa)
axis tight
legend('error for NTPFAlin','error for NTPFAopt');
hold off

% figure(4)
% clf
% x = G.cells.centroids(:, 1);
% plot(x,p_exact,'.')
% hold on
% plot(x,state2.pressure,'+');
% plot(x,state.pressure,'d');
% %plot(x,mpfa.pressure,'.');
% %legend('ntpfaOPT','ntpfaLIN','MPFA');
% legend('exact solution','NTPFA','NTPFAopt');

