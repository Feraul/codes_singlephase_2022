function [ netas,coeff] =netas_Interp_LPEW(O, P, T, Qo, ni,kmap,gravaux)

global esurn2 elem
% Retorna os netas.
% esta variavel podemos achar no pag. 07 do artigo  chinês

%Prealocação do vetore.%
netas=zeros(esurn2(ni+1)-esurn2(ni),2);
R=[0 1 0; -1 0 0; 0 0 0];
%Loop que percorre os elementos em torno do nó "ni".%
K=zeros(3);
for k=1:size(netas,1),
    
    %Essa é UMA maneira de construir os tensores
    K(1,1) = kmap(elem(k,5),2);
    K(1,2) = kmap(elem(k,5),3);
    K(2,1) = kmap(elem(k,5),4);
    K(2,2) = kmap(elem(k,5),5);
    %Preenchimento da segunda coluna do vetor "netas".%
    
    if (k==size(netas,1))&&(size(P,1)==size(O,1))
        v1=O(k,:)-Qo;
        v2=P(1,:)-Qo;
        ce=cross(v1,v2); % produto vetorial
        H2=norm(ce)/norm(v2); % altura
        netas(k,2)=norm(T(1,:)-Qo)/H2;
        coeff(k,2)=gravaux(k,2);
    else
        v1=O(k,:)-Qo;
        v2=P(k+1,:)-Qo;
        ce=cross(v1,v2);
        H2=norm(ce)/norm(v2); % altura
        netas(k,2)=norm(T(k+1,:)-Qo)/H2;
        coeff(k,2)=gravaux(k,2);
    end
    
    %%%%%%%Fim do Preenchimento da segunda coluna do vetor "netas".%%%%%%%
    
    
    %Preenchimento da primeira coluna do vetor "netas".%
    
    v1=O(k,:)-Qo;
    v2=P(k,:)-Qo;
    ce=cross(v1,v2); % produto vetorial
    H1=norm(ce)/norm(v2); % altura
    netas(k,1)=norm(T(k,:)-Qo)/H1;
    coeff(k,1)=gravaux(k,1);
    %%%%%Fim do Preenchimento da primeira coluna do vetor "netas".%%%%%%%

end

