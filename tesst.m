 G = computeGeometry(cartGrid([10,10], [1,1]));
G = computeGeometry(twister(G));
srci = @(x) -sin(x);

src_ex = zeros(G.cells.num,1);

%p_ex = @(cntr) sin(pi*cntr(1))*sin(pi*(cntr(2)+1/2))+1;
p_ex = @(cntr) sin(cntr(1)) + 2*cntr(2);
p_exact = zeros(G.cells.num,1);
for j=1:G.cells.num
    p_exact(j) = p_ex(G.cells.centroids(j,:));
    src_ex(j) = srci(G.cells.centroids(j,1));
end

rock = makeRock(G, 100*milli*darcy, 0.002);
K = [1, 0.5; 0.5, 1];
rock.perm = repmat([K(1,1), K(1,2), K(2,2)], G.cells.num, 1);
fluid = initSimpleADIFluid();

[src,W,bc] = deal([]);
%src = addSource(src, [1:G.cells.num],src_ex,'sat',[1]);
bndryFaces=boundaryFaces(G);
cntBndryFaces = G.faces.centroids(bndryFaces,:);
for k = 1 : length(bndryFaces)
    p_ex_temp = p_ex(cntBndryFaces(k,:));
    bc = addBC(bc, bndryFaces(k), 'pressure', p_ex_temp,'sat',[1]);
end

state0 = initResSol(G, 0, [0 1]);
model = PressureOilWaterModelNTPFAopt(G,rock,fluid);
state = incompSinglePhaseNTPFA(model, state0,'bc', bc, 'src',src);

error = sqrt((sum((G.cells.volumes.*(p_exact-state.pressure)).^2))/(sum((G.cells.volumes.*p_exact).^2)));


%profile on / off / viewer