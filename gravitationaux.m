function [G,g]=gravitationaux(kmap,gravelem)
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
    ve1=0.5*(coord(bedge(ifacont,2),:)-coord(bedge(ifacont,1),:));

    % tensor de permeabilidade do elemento a esquerda
    Klef(1,1)=kmap(elem(lef,5),2);
    Klef(1,2)=kmap(elem(lef,5),3);
    Klef(2,1)=kmap(elem(lef,5),4);
    Klef(2,2)=kmap(elem(lef,5),5);

    Keq=Klef;

    g(ifacont,1)=dot((R*ve1')'*Keq,(gravelem(lef,:)));

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
    vd1=0.5*(coord(inedge(iface,2),1:2)-coord(inedge(iface,1),1:2));
    Keq=inv((dj1*inv(K3)+dj2*inv(K4))); % equation 21
    graveq=((dj1*gravelem(lef,1:2)+dj2*gravelem(rel,1:2))'); % equation 22
    g(iface+size(bedge,1),1)=dot(((R1*vd1')')*Keq, graveq);% equation 20

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
