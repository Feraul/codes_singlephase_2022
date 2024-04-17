function [G,g]=gravitation(kmap,gradgravelem,gravface,Hesq, Kde,Kn,Kt,Ded,...
    grav_elem_escalar,gravno,w,nflagno,wg)
global inedge bedge elem centelem coord normals strategy esurn1 esurn2 benchmark
Klef=zeros(3,3);
G=zeros(size(elem,1),1);
R=[0 1 0;-1 0 0 ; 0 0 0];
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

    O=centelem(lef,:); % baricentro do elemento a esuqerda
    B1=bedge(ifacont,1);
    B2=bedge(ifacont,2);
    nor=norm(coord(B1,:)-coord(B2,:));
    if  strcmp(strategy,'starnoni')

        g(ifacont,1)=dot((R*ve1')'*Klef,(gradgravelem(lef,:)));

    elseif strcmp(strategy,'GravConsist')
        if bedge(ifacont,5)>200
            g(ifacont,1)=0;
        else
            if strcmp(benchmark,'starnonigrav1')
                 A=(Kn(ifacont)/(Hesq(ifacont)*nor));
                gravno1=gravno(bedge(ifacont,1),1);
                gravno2=gravno(bedge(ifacont,2),1);
%         g(ifacont,1)=-A*(((O-coord(B2,:)))*(coord(B1,:)-coord(B2,:))'*gravno1+...
%             (O-coord(B1,:))*(coord(B2,:)-coord(B1,:))'*gravno2-(nor^2)*grav_elem_escalar(lef))...
%             -(gravno2-gravno1)*Kt(ifacont);

                g(ifacont,1)=-dot((R*ve1')'*Klef,(gradgravelem(lef,:)));
                
            else
                g(ifacont,1)=-dot((R*ve1')'*Klef,(gradgravelem(lef,:)));
            end
        end

    else
        %g(ifacont,1)=dot((R*ve1')'*Klef,(gravface(ifacont,:)));
        g(ifacont,1)=dot((R*ve1')'*Klef,(gradgravelem(lef,:)));
    end
    G(lef,1)=G(lef,1)-g(ifacont,1);

end
for iface=1:size(inedge,1)
    % elementos a esquerda e a direita
    lef=inedge(iface,3);
    rel=inedge(iface,4);

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
        vd1=coord(inedge(iface,2),1:3)-coord(inedge(iface,1),1:3);
        Keq=inv((dj1*inv(K3)+dj2*inv(K4))); % equation 21
        graveq=((dj1*gradgravelem(lef,1:3)+dj2*gradgravelem(rel,1:3))'); % equation 22
        g(iface+size(bedge,1),1)=dot(((R1*vd1')')*Keq, graveq);% equation 20
    elseif strcmp(strategy,'GravConsist')
        %Determinação dos centróides dos elementos à direita e à esquerda.
        %====================================================================================
        no1=inedge(iface,1);
        no2=inedge(iface,2);
        C1 = centelem(inedge(iface,3),1:3); % baricentro do elemento a esquerda
        C2 = centelem(inedge(iface,4),1:3); % baricentro do elemento direito
        vcen = C2 - C1;
        vd1 = coord(inedge(iface,2),1:3) - coord(inedge(iface,1),1:3);
        ve2 = C1 - coord(inedge(iface,1),1:3);
        vd2 = C2 - coord(inedge(iface,1),1:3);     %Do início da aresta até o
        %centro da célula da direita.

        ce = cross(vd1,ve2);
        H1 = norm(ce)/norm(vd1); % altura a esquerda
        %Determinação das alturas dos centróides dos elementos à direita e à%
        %esquerda.                                                          %
        cd = cross(vd1,vd2);
        H2 = norm(cd)/norm(vd1); % altura a direita
        % tensor de permeabilidade do elemento a esquerda

        K3(1,1)=kmap(elem(lef,5),2);
        K3(1,2)=kmap(elem(lef,5),3);
        K3(2,1)=kmap(elem(lef,5),4);
        K3(2,2)=kmap(elem(lef,5),5);

        % tensor de permeabilidade do elemento a direita

        K4(1,1)=kmap(elem(rel,5),2);
        K4(1,2)=kmap(elem(rel,5),3);
        K4(2,1)=kmap(elem(rel,5),4);
        K4(2,2)=kmap(elem(rel,5),5);
        vd11=coord(inedge(iface,2),1:2) - coord(inedge(iface,1),1:2);
        % calculo das constantes tangenciais e normais em cada face interna
        Kn1 = (RotH(vd11)'*K3*RotH(vd11))/norm(vd11)^2;
        Kt1 = (RotH(vd11)'*K3*(vd11)')/norm(vd11)^2;

        Kn2 = (RotH(vd11)'*K4*RotH(vd11))/norm(vd11)^2;
        Kt2 = (RotH(vd11)'*K4*(vd11)')/norm(vd11)^2;

        Kde1 = -((Kn1*Kn2))/(Kn1*H2 + Kn2*H1);
        % Ded: constante que tem constantes geometricas + contantes
        % tangeciais
        Ded1 = (dot(vd1,vcen)/norm(vd1)^2) -...
            (1/norm(vd1))*((Kt2/Kn2)*H2 + (Kt1/Kn1)*H1);
        %aproximacao do termo grav do elemento L na face IJ
        gleft=dot(RotH(vd11)'*K3,(gradgravelem(lef,1:2)));
        %aproximacao do termo grav do elemento R na face JI
        gright=dot(-RotH(vd11)'*K4,(gradgravelem(rel,1:2)));

        %if nflagno(no1,1)>200
           g1=wg(no1);
        %else
        %   g1=gravno(no1,1);
        %end

        %if nflagno(no2,1)>200
           g2=wg(no2);
        %else
        %   g2=gravno(no2,1);
        %end
%         g1=0;
%         nec1=esurn2(no1+1)- esurn2(no1);
%         nec2=esurn2(no2+1)- esurn2(no2);
%         %if nflagno(no1,1)>200
%             for j=1:nec1
%                 element1=esurn1(esurn2(no1)+j);
%                 %if element1~=lef && element1~=rel
%                     g1=g1+wg(esurn2(no1)+j);
%                 %end
%             end

        %else
        %    g1=gravno(no1,1);
        %end
        %g2=0;

        %if nflagno(no2,1)>200
        %    for j=1:nec2
        %        element2=esurn1(esurn2(no2)+j);
        %        %if element2~=lef && element2~=rel
        %            g2=g2+wg(esurn2(no2)+j);
        %        %end
        %    end
        %else
        %    g2=gravno(no2,1);
        %end

        vd1=coord(inedge(iface,2),1:2)-coord(inedge(iface,1),1:2);
        Keq=inv((dj1*inv(K3)+dj2*inv(K4))); % equation 21
        graveq=((dj1*gradgravelem(lef,1:2)+dj2*gradgravelem(rel,1:2))'); % equation 22

        %g(iface+size(bedge,1),1)=-dot((RotH(vd11)')*Keq, graveq)+Kde1*Ded1*norm(vd1)*(g1-g2);
        g(iface+size(bedge,1),1)=-Kde1*((H2/Kn2)*gright-(H1/Kn1)*gleft)+Kde1*Ded1*norm(vd1)*(g1-g2);
        %g(iface+size(bedge,1),1)=Kde1*norm(vd1)*(grav_elem_escalar(rel)-grav_elem_escalar(lef))+Kde1*Ded1*norm(vd1)*(g1-g2);
    end
    G(lef,1)=G(lef,1)-g(iface+size(bedge,1),1);
    G(rel,1)=G(rel,1)+g(iface+size(bedge,1),1);
end
end

function [RH]=RotH(vi)
%Função que retorna a rotação de um certo vetor em 90 graus,
%horariamente, na primeira coluna da matriz R, e horariamente,
%na segunda. Restringe um vetor de 3 coordenadas a um de 2, considerando
%que a terceira coordenada é nula.
% vi2=zeros(2,1);
% if size(vi)~=[3 1]
%     vi=vi';
%     vi2(1)=vi(1);
%     vi2(2)=vi(2);
% end
RH=[0 1;-1 0 ]*vi';
end
