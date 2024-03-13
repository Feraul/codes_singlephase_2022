function [G,g]=gravitationHP(kmap,grav_elem_escalar,gravface,...
    pointarmonic,weightDMP,gravelem,parameter,nflagface)
global inedge bedge elem centelem coord strategy benchmark
Klef=zeros(3,3);
G=zeros(size(elem,1),1);
R=[0 1 0;-1 0 0 ; 0 0 0];
K1=zeros(3,3);
K2=zeros(3,3);
K3=zeros(2,2);
K4=zeros(2,2);
[gravaux]=gravinterpHP(weightDMP,parameter,gravface,grav_elem_escalar,nflagface);

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
    norma=norm(ve1);
    if  strcmp(strategy,'starnoni')
%         if bedge(ifacont,5)>200
%             g(ifacont,1)=dot((R*ve1')'*Klef,(gravface(ifacont,:)));
%         else
            
            g(ifacont,1)=dot((R*ve1')'*Klef,(gravelem(lef,:)));
       %end
    elseif strcmp(strategy,'GravConsist')
        if bedge(ifacont,5)>200
            g(ifacont,1)=0;
        else
            if strcmp(benchmark,'starnonigrav1')
                % g(ifacont,1)=dot((R*ve1')'*Klef,(gravelem(lef,:)));
                % calculo dos fluxo parcial a esquerda
                % os nós que conforman os pontos de interpolação no elemento a esquerda
                auxfacelef1=parameter(1,3,ifacont);
                auxfacelef2=parameter(1,4,ifacont);
                gaux1=gravaux(auxfacelef1);
                gaux2=gravaux(auxfacelef2);
                g(ifacont,1)=norma*((parameter(1,1,ifacont)+parameter(1,2,ifacont))*grav_elem_escalar(lef)-...
                    parameter(1,1,ifacont)*gaux1-parameter(1,2,ifacont)*gaux2);
            else
                g(ifacont,1)=-dot((R*ve1')'*Klef,(gravelem(lef,:)));
            end
        end
    end
    G(lef,1)=G(lef,1)-g(ifacont,1);
    
end
for iface=1:size(inedge,1)
    % elementos a esquerda e a direita
    lef=inedge(iface,3);
    rel=inedge(iface,4);
    if  strcmp(strategy,'starnoni')
        % calculo do ponto meio da face
        vm=(coord(inedge(iface,1),:)+coord(inedge(iface,2),:))*0.5;
        
        % calculo da distancia do centro ao ponto meio da face
        dj1=norm(centelem(lef,:)-vm);
        dj2=norm(centelem(rel,:)-vm);
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
        
        %=======================================================================
        
        vd1=coord(inedge(iface,2),:)-coord(inedge(iface,1),:);
        norma=norm(vd1);
        ifactual=iface+size(bedge,1);
        % os nós que conforman os pontos de interpolação no elemento a esquerda
        auxfacelef1=parameter(1,3,ifactual);
        auxfacelef2=parameter(1,4,ifactual);
        % os nós que conforman os pontos de interpolação no elemento a direita
        auxfacerel1=parameter(2,3,ifactual);
        auxfacerel2=parameter(2,4,ifactual);
        % calculo dos fluxo parcial a esquerda
        gaux1=gravaux(auxfacelef1);
        gaux2=gravaux(auxfacelef2);
        % fluxo gravitacional: esquerda
        fluxesq=norma*((parameter(1,1,ifactual)+parameter(1,2,ifactual))*grav_elem_escalar(lef)-...
            parameter(1,1,ifactual)*gaux1-parameter(1,2,ifactual)*gaux2);
        % calculo dos fluxo parcial direita
        
        gaux3=gravaux(auxfacerel1);
        gaux4=gravaux(auxfacerel2);
        % fluxo gravitacional: direita
        fluxdireit=norma*((parameter(2,1,ifactual)+parameter(2,2,ifactual))*grav_elem_escalar(rel)-...
            parameter(2,1,ifactual)*gaux3-parameter(2,2,ifactual)*gaux4);
        % Calculo das contribuições do elemento a esquerda
        mulef=weightDMP(ifactual-size(bedge,1),1);
        murel=weightDMP(ifactual-size(bedge,1),2);
        % calculo do fluxo unico na face
        g(iface+size(bedge,1),1)= (murel*fluxesq-mulef*fluxdireit);
    end
    G(lef,1)=G(lef,1)-g(iface+size(bedge,1),1);
    G(rel,1)=G(rel,1)+g(iface+size(bedge,1),1);
end
end
