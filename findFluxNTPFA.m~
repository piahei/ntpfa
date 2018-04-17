function  v = findFluxNTPFA(model,p)

op = model.operators;

L_ij1 = op.L_ij1;
L_ij2 = op.L_ij2;
L_ji1 = op.L_ji1;
L_ji2 = op.L_ji2;
A_ij = op.A_ij;
A_ji = op.A_ji;
pos_ij = op.pos_ij;
pos_ji = op.pos_ji;

v = initVariablesADI(zeros(numel(p.val),1));
v = op.Grad(v);

for iFace=1:numel(op.intFace)
    ci = model.G.faces.neighbors(op.intFace(iFace),1); %neighbor i
    Fi = model.G.cells.faces(model.G.cells.facePos(ci):model.G.cells.facePos(ci+1)-1,1); %faces i
    cj = model.G.faces.neighbors(op.intFace(iFace),2); %neighbor j
    Fj = model.G.cells.faces(model.G.cells.facePos(cj):model.G.cells.facePos(cj+1)-1,1); %faces j
    
    % Remove common face and external face:
    ni = model.G.faces.neighbors(Fi,:); ni(ni==ci) = []; %neighbors to i
    nj = model.G.faces.neighbors(Fj,:); nj(nj==cj) = []; %neighbors to j
    zi = find(~ni); zj = find(~nj); %external faces
   
    
    %remove neighbor with common face iFace
    pos_ci = find(ni==cj);  %find position of cj in ni
    pos_cj = find(nj==ci);  %find position of ci in nj
    ni=nonzeros(ni); nj=nonzeros(nj); %remove external "cell"
    nj(nj==ci)=[];
    ni(ni==cj)=[];
    
    %positions in lambda and positions to be removed from lambda
    ij_temp = pos_ij(iFace):pos_ij(iFace+1)-1;
    ij_temp([pos_ci,zi])=[];  
    ji_temp = pos_ji(iFace):pos_ji(iFace+1)-1;
    ji_temp([pos_cj,zj])=[];
    
    L_ij = sum(L_ij1(ij_temp))*p.val(ci)+sum(L_ij2(ij_temp).*p.val(ni));
    L_ji = sum(L_ji1(ji_temp))*p.val(cj)+sum(L_ji2(ji_temp).*p.val(nj));
    
    if L_ij == 0 && L_ji == 0
        mu_ij = 0.5;
        mu_ji = 0.5;
    else
        mu_ij = abs(L_ji)./(abs(L_ij)+abs(L_ji));
        mu_ji = abs(L_ij)./(abs(L_ij)+abs(L_ji));
    end
    
    % Compute flux
    % forskjellig r? Dirichlet?
    
    r_ij = mu_ji*L_ji-mu_ij*L_ij; %zero if L_ij*L_ji>=0
    
    ijPos = pos_ij(iFace):pos_ij(iFace+1)-1;
    jiPos = pos_ji(iFace):pos_ji(iFace+1)-1;
    
    %transmisibilities
    T_ij = mu_ij*sum(A_ij(ijPos))+mu_ji*A_ji(jiPos(pos_cj));
    T_ji = mu_ji*sum(A_ji(jiPos))+mu_ij*A_ij(ijPos(pos_ci));
    
    
    %flux
    v(iFace) = (T_ij+(abs(r_ij)+r_ij)./(2*(p(ci)+eps)))*p(ci)-...
        (T_ji+(abs(r_ij)-r_ij)./(2*(p(cj)+eps)))*p(cj);     
    
end
   
   
end