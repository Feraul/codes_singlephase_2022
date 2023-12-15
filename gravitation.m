function [G,g]=gravitation(kmap,gravelem,gravface,Hesq, Kde,Kn,Kt,Ded,...
    grav_elem_escalar,gravno,w,nflagno,wg)
global inedge bedge elem centelem coord normals strategy esurn1 esurn2
Klef=zeros(3,3);
G=zeros(size(elem,1),1);
R=[0 1 0;-1 0 0 ; 0 0 0];
K1=zeros(3,3);
K2=zeros(3,3);
K3=zeros(2,2);
K4=zeros(2,2);
for ifacont=1:size(bedge,1)
    % elemento a esquerda
    lef=bedge(ifacont,3);
    % vetor de orientacao da face em questao
    ve1=coord(bedge(ifacont,2),:)-coord(bedge(ifacont,1),:);
    
    % tensor de permeabilidade do elemento a esquerda
    Klef(1,1)=kmap(elem(lef,5),2);
    Klef(1,2)=kmap(elem(lef,5),3);
    Klef(2,1)=kmap(elem(lef,5),4);
    Klef(2,2)=kmap(elem(lef,5),5);
    
    Keq=Klef;
    if  strcmp(strategy,'starnoni')
        if bedge(ifacont,5)>200
            g(ifacont,1)=dot((R*ve1')'*Keq,(gravface(ifacont,:)));
        else
            g(ifacont,1)=dot((R*ve1')'*Keq,(gravelem(lef,:)));
        end
    elseif strcmp(strategy,'GravConsist')
        g(ifacont,1)=dot((R*ve1')'*Keq,(gravelem(lef,:)));
    elseif strcmp(strategy,'inhouse')
        v0=coord(bedge(ifacont,2),:)-coord(bedge(ifacont,1),:); %fase.
        v1=centelem(bedge(ifacont,3),:)-coord(bedge(ifacont,1),:);
        v2=centelem(bedge(ifacont,3),:)-coord(bedge(ifacont,2),:);
        normcont=norm(v0);
        
        % Tratamento do nó nos vértices 2 e 4%
        A=-Kn(ifacont)/(Hesq(ifacont)*norm(v0));
        
        if bedge(ifacont,5)<200
            g1=gravno(bedge(ifacont,1),1); % gravidade no vertice 1
            g2=gravno(bedge(ifacont,2),1); % gravidade no vertice 2
        else
            no1=bedge(ifacont,1);
            no2=bedge(ifacont,2);
            
            nec1=esurn2(no1+1)- esurn2(no1);
            nec2=esurn2(no2+1)- esurn2(no2);
            g1=0;
            if nflagno(no1,1)>200
                for j=1:nec1
                    element1=esurn1(esurn2(no1)+j);
                    g1=g1+w(esurn2(no1)+j)*grav_elem_escalar(element1);
                end
            else
                g1=gravno(no1,1);
            end
            g2=0;
            if nflagno(no2,1)>200
                for j=1:nec2
                    element2=esurn1(esurn2(no2)+j);
                    g2=g2+w(esurn2(no2)+j)*grav_elem_escalar(element2);
                end
            else
                g2=gravno(no2,1);
            end
        end
        g(ifacont,1)=-(A*(dot(v2,-v0)*g1+dot(v1,v0)*g2-(normcont^2*grav_elem_escalar(lef)))-(g2-g1)*Kt(ifacont));
        
    else
        %g(ifacont,1)=dot((R*ve1')'*Keq,(gravface(ifacont,:)));
        g(ifacont,1)=dot((R*ve1')'*Keq,(gravelem(lef,:)));
    end
    G(lef,1)=G(lef,1)-g(ifacont,1);
    
end
for iface=1:size(inedge,1)
    % elementos a esquerda e a direita
    lef=inedge(iface,3);
    rel=inedge(iface,4);
    % vertices
    no1=inedge(iface,1);
    no2=inedge(iface,2);
    
    % calculo do ponto meio da face
    vm=(coord(inedge(iface,1),:)+coord(inedge(iface,2),:))*0.5;
    
    % calculo da distancia do centro ao ponto meio da face
    dj1=norm(centelem(lef,:)-vm);
    dj2=norm(centelem(rel,:)-vm);
    
    
    if  strcmp(strategy,'starnoni')
        % tensor de permeabilidade do elemento a esquerda
        R1=[0 1 ;-1 0 ];
        K3(1,1)=kmap(elem(lef,5),2);
        K3(1,2)=kmap(elem(lef,5),3);
        K3(2,1)=kmap(elem(lef,5),4);
        K3(2,2)=kmap(elem(lef,5),5);
        
        % tensor de permeabilidade do elemento a direita
        K4(1,1)=kmap(elem(rel,5),2);
        K4(1,2)=kmap(elem(rel,5),3);
        K4(2,1)=kmap(elem(rel,5),4);
        K4(2,2)=kmap(elem(rel,5),5);
        vd1=coord(inedge(iface,2),1:2)-coord(inedge(iface,1),1:2);
        Keq=inv((dj1*inv(K3)+dj2*inv(K4))); % equation 21
        graveq=((dj1*gravelem(lef,1:2)+dj2*gravelem(rel,1:2))'); % equation 22
        g(iface+size(bedge,1),1)=dot(((R1*vd1')')*Keq, graveq);% equation 20
    elseif strcmp(strategy,'GravConsist')
        %Determinação dos centróides dos elementos à direita e à esquerda.%
        C1 = centelem(inedge(iface,3),:); % baricentro do elemento a esquerda
        C2 = centelem(inedge(iface,4),:); % baricentro do elemento direito
        vcen = C2 - C1;
        vd1 = coord(inedge(iface,2),:) - coord(inedge(iface,1),:);
        ve2 = C1 - coord(inedge(iface,1),:);
        vd2 = C2 - coord(inedge(iface,1),:);     %Do início da aresta até o
        %centro da célula da direita.
        
        ce = cross(vd1,ve2);
        H1 = norm(ce)/norm(vd1); % altura a esquerda
        %Determinação das alturas dos centróides dos elementos à direita e à%
        %esquerda.                                                          %
        cd = cross(vd1,vd2);
        H2 = norm(cd)/norm(vd1); % altura a direita
        % tensor de permeabilidade do elemento a esquerda
        
        K1(1,1)=kmap(elem(lef,5),2);
        K1(1,2)=kmap(elem(lef,5),3);
        K1(2,1)=kmap(elem(lef,5),4);
        K1(2,2)=kmap(elem(lef,5),5);
        
        % tensor de permeabilidade do elemento a direita
        
        K2(1,1)=kmap(elem(rel,5),2);
        K2(1,2)=kmap(elem(rel,5),3);
        K2(2,1)=kmap(elem(rel,5),4);
        K2(2,2)=kmap(elem(rel,5),5);
        % calculo das constantes tangenciais e normais em cada face interna
        Kn1 = (RotH(vd1)'*K1*RotH(vd1))/norm(vd1)^2;
        Kt1 = (RotH(vd1)'*K1*(vd1)')/norm(vd1)^2;
        
        Kn2 = (RotH(vd1)'*K2*RotH(vd1))/norm(vd1)^2;
        Kt2 = (RotH(vd1)'*K2*(vd1)')/norm(vd1)^2;
        
        Kde = ((Kn1*Kn2))/(Kn1*H2 + Kn2*H1);
        % Ded: é uma constante que tem constantes geometricas + contantes
        % tangeciais
        Ded = (dot(vd1,vcen)/norm(vd1)^2) -(1/norm(vd1))*((Kt2/Kn2)*H1 + (Kt1/Kn1)*H2);
        
        gleft=dot(RotH(vd1)'*K1,(gravelem(lef,:)));
        gright=dot(RotH(-vd1)'*K2,(gravelem(rel,:)));
        
        if nflagno(no1,1)>200
            g1=wg(no1);
        else
            g1=gravno(no1,1);
        end
        
        if nflagno(no2,1)>200
            g2=wg(no2);
        else
            g2=gravno(no2,1);
        end
        %------------------------------------------------------------------
        g(iface+size(bedge,1),1)=-(Kde*((H2/Kn2)*gright-(H1/Kn1)*gleft-Ded*(g2-g1)));
    else
        
        no1=inedge(iface,1);
        no2=inedge(iface,2);
        g1=0;
        nec1=esurn2(no1+1)- esurn2(no1);
        nec2=esurn2(no2+1)- esurn2(no2);
        if nflagno(no1,1)>200
            for j=1:nec1
                element1=esurn1(esurn2(no1)+j);
                g1=g1+w(esurn2(no1)+j)*grav_elem_escalar(element1);
            end
        else
            g1=gravno(no1,1);
        end
        g2=0;
        
        if nflagno(no2,1)>200
            for j=1:nec2
                element2=esurn1(esurn2(no2)+j);
                g2=g2+w(esurn2(no2)+j)*grav_elem_escalar(element2);
            end
        else
            g2=gravno(no2,1);
        end
        g(iface+size(bedge,1),1)= -Kde(iface)*(grav_elem_escalar(rel)-grav_elem_escalar(lef)-Ded(iface)*(g2-g1));
    end
    G(lef,1)=G(lef,1)-g(iface+size(bedge,1),1);
    G(rel,1)=G(rel,1)+g(iface+size(bedge,1),1);
end
end

function [RH]=RotH(vi)
%Função que retorna a rotação de um certo vetor em 90 graus,
%anti-horariamente, na primeira coluna da matriz R, e horariamente,
%na segunda. Restringe um vetor de 3 coordenadas a um de 2, considerando
%que a terceira coordenada é nula.
% vi2=zeros(2,1);
% if size(vi)~=[3 1]
%     vi=vi';
%     vi2(1)=vi(1);
%     vi2(2)=vi(2);
% end
RH=[0 1 0;-1 0 0;0 0 0]*vi';
end
