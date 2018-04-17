function [L_ij1,L_ij2,L_ji1,L_ji2,A_ij,A_ji,pos_ij,pos_ji,intFace] = ntpfa(G,rock)

noFaces = G.faces.num;   %number of faces;
numFaceCell = diff(G.cells.facePos); %number of faces for each cell
%faceAreas = G.faces.areas;  

intFace = zeros(noFaces,1);  
for l=1:noFaces
    intFaceNo_temp = G.faces.neighbors(l,:);
    if intFaceNo_temp(1)==0 || intFaceNo_temp(2)==0
        intFace(l) = 0;
    else
        intFace(l) = l;
    end
end
intFace = nonzeros(intFace);  %delete edge faces
numIntFace = numel(intFace);

%Find harmonic averaging point (hap) for each face
hap = zeros(G.faces.num,2);
 
beta_ij = zeros(G.faces.num,1); %beta for internal faces (zero for external)
beta_ji = zeros(G.faces.num,1);
O_i = zeros(G.faces.num,1); %omega_ij  
O_j = zeros(G.faces.num,1); %omega_ji
ii=1;
for i = 1:noFaces  %over all faces
    n1 = G.faces.neighbors(i,1); %neighbor cell 1
    n2 = G.faces.neighbors(i,2); %neighbor cell 2
    node1 = G.faces.nodes(ii);
    node2 = G.faces.nodes(ii+1);
    fn1 = G.nodes.coords(node1,:);  %coordinates to node 1
    fn2 = G.nodes.coords(node2,:);  %coordinates to node 2
    ii = ii+2;
    
    if n1 == 0 || n2 == 0  %external face -> hap = face centroid
        hap(i,:) = G.faces.centroids(i,:);
     
    else  %internal faces
        K1 = rock.perm(n1);
        K2 = rock.perm(n2);
        centroid1 = G.cells.centroids(n1,:); %centroid of neighboring cell 1
        centroid2 = G.cells.centroids(n2,:); %centroid of neighboring cell 2
        dist_i = abs((fn2(2)-fn1(2))*centroid1(1)-(fn2(1)-fn1(1))*centroid1(2)+fn2(1)*fn1(2)-fn2(2)*fn1(1))/...
        sqrt((fn2(2)-fn1(2))^2+(fn2(1)-fn1(1))^2);    % shortest distance from cell centroid 1 to face i
        dist_j = abs((fn2(2)-fn1(2))*centroid2(1)-(fn2(1)-fn1(1))*centroid2(2)+fn2(1)*fn1(2)-fn2(2)*fn1(1))/...
        sqrt((fn2(2)-fn1(2))^2+(fn2(1)-fn1(1))^2);    % shortest distance from cell centroid 2 to face i
        beta_ij(i) = K1*(G.faces.normals(i,:)*G.faces.normals(i,:)')/dist_i;
        %faceAreas(i)^2*
        beta_ji(i) = K2*(G.faces.normals(i,:)*G.faces.normals(i,:)')/dist_j;
        %faceAreas(i)^2*
        O_i(i) = beta_ij(i)/(beta_ij(i)+beta_ji(i));  %omega
        O_j(i) = beta_ji(i)/(beta_ij(i)+beta_ji(i));
        hap(i,:) = (beta_ij(i)+beta_ji(i))\...
          (beta_ij(i)*centroid1'+beta_ji(i)*centroid2'+(K1-K2)*G.faces.normals(i,:)');
    
    end
    
end

%Preallocate vectors for alpha*omega and lambda, 
%Split lambda in two

L_ij1 = zeros(max(numFaceCell)*numIntFace,1); %A_ij*omega_ik  (Lambda)
L_ji1 = zeros(max(numFaceCell)*numIntFace,1); %A_ji*omega_jk  (Lambda)
L_ij2 = zeros(max(numFaceCell)*numIntFace,1); %A_ij*omega_ki
L_ji2 = zeros(max(numFaceCell)*numIntFace,1); %A_ji*omega_kj 
 
%alpha_ij = zeros(max(numFaceCell)*numIntFace,1);
A_ij = zeros(max(numFaceCell)*numIntFace,1); %alpha_ij,k*omega_ki
%alpha_ji = zeros(max(numFaceCell)*numIntFace,1);
A_ji = zeros(max(numFaceCell)*numIntFace,1); %alpha_ji,k*omega_kj
pos_ij = ones(numIntFace+1,1);
pos_ji = ones(numIntFace+1,1);
pi=1; pj=1; 
%Find decomposed vector using hap for internal faces
  for iFace=1:numIntFace %number of internal faces
      cell_i = G.faces.neighbors(intFace(iFace),1); %neighbor 1
      cell_j = G.faces.neighbors(intFace(iFace),2); %neighbor 2
      Fi = G.cells.faces(G.cells.facePos(cell_i):G.cells.facePos(cell_i+1)-1,1); %faces to cell i
      Fj = G.cells.faces(G.cells.facePos(cell_j):G.cells.facePos(cell_j+1)-1,1); %faces to cell j
      faceNormal_ij = rock.perm(cell_i)*G.faces.normals(intFace(iFace),:); %normal from cell 1->2
      faceNormal_ji = -rock.perm(cell_j)*G.faces.normals(intFace(iFace),:); %normal from cell 2->1
      
      pos_ij(iFace+1) = pi+length(Fi); pi = pi+length(Fi);
      pos_ji(iFace+1) = pj+length(Fj); pj = pj+length(Fj);
      
      %find alpha's
      [alpha_ij,~,~] = decomp(faceNormal_ij,Fi,G.cells.centroids(cell_i,:),hap(Fi,:));
      [alpha_ji,~,~] = decomp(faceNormal_ji,Fj,G.cells.centroids(cell_j,:),hap(Fj,:)); 
       
      A_ij_temp = alpha_ij.*O_j(Fi);
      A_ji_temp = alpha_ji.*O_j(Fj);
      
      L_ij1(pos_ij(iFace):pos_ij(iFace+1)-1) = A_ij_temp.*O_i(Fi);
      L_ij2(pos_ij(iFace):pos_ij(iFace+1)-1) = A_ij_temp.*O_j(Fi);
      L_ji1(pos_ji(iFace):pos_ji(iFace+1)-1) = A_ji_temp.*O_i(Fj);
      L_ji2(pos_ji(iFace):pos_ji(iFace+1)-1) = A_ji_temp.*O_j(Fj); 
      
      A_ij(pos_ij(iFace):pos_ij(iFace+1)-1) = A_ij_temp;
      A_ji(pos_ij(iFace):pos_ij(iFace+1)-1) = A_ji_temp;
  end
  L_ij1 = L_ij1(1:pos_ij(end)-1); L_ij2 = L_ij2(1:pos_ij(end)-1);  %remove zero elements
  L_ji1 = L_ji1(1:pos_ji(end)-1); L_ji2 = L_ji2(1:pos_ji(end)-1);  %remove zero elements
  A_ij = A_ij(1:pos_ij(end)-1);
  A_ji = A_ji(1:pos_ji(end)-1);
end

function [alpha,t,d] = decomp(faceNormal,facesCell,centroid,hap)

noFace = length(facesCell);
alpha = zeros(noFace,1);
t = bsxfun(@minus,centroid,hap)'; %vec with x_i-x_k

funk = @(x) sum(x);
option = optimset('Display','off');
alpha = fmincon(funk,alpha,[],[],t,faceNormal',zeros(noFace,1),Inf*ones(noFace,1),[],option);

d = t*alpha;

end