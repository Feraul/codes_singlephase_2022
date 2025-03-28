function [ netas,gaux] =netas_Interp_LPEW(O, P, T, Qo, ni,kmap,gravelem,flagbedge)

global esurn2 elem esurn1 gravitational
% Retorna os netas.
% esta variavel podemos achar no pag. 07 do artigo  chinês

%Prealocação do vetore.%
netas=zeros(esurn2(ni+1)-esurn2(ni),2);
gaux=zeros(esurn2(ni+1)-esurn2(ni),6);

% matriz de rotacao em sentido horario
R=[0 1 0; -1 0 0; 0 0 0];
%Loop que percorre os elementos em torno do nó "ni".%
K=zeros(3);

for k=1:size(netas,1)
    j=esurn1(esurn2(ni)+k);
    %Essa é UMA maneira de construir os tensores
    K(1,1) = kmap(elem(j,5),2);
    K(1,2) = kmap(elem(j,5),3);
    K(2,1) = kmap(elem(j,5),4);
    K(2,2) = kmap(elem(j,5),5);
    %Preenchimento da segunda coluna do vetor "netas".%

    if (k==size(netas,1))&&(size(P,1)==size(O,1))
        v1=O(k,:)-Qo;
        v2=P(1,:)-Qo;
        ce=cross(v1,v2); % produto vetorial
        H2=norm(ce)/norm(v2); % altura
        netas(k,2)=norm(T(1,:)-Qo)/H2;
        if strcmp(gravitational,'yes')
            if flagbedge(k,1)~=0
                % para vertice no contorno do dominio
                if flagbedge(k,2)>200
                    gaux(k,4:6)=(R*(Qo-T(1,:))')';
                else
                    gaux(k,4:6)=(R*(Qo-T(1,:))')';
                end
            else
                % para vertice no interior do dominio
                gaux(k,4:6)=(R*(Qo-T(1,:))')';
            end
        end
    else
        v1=O(k,:)-Qo;
        v2=P(k+1,:)-Qo;
        ce=cross(v1,v2);
        H2=norm(ce)/norm(v2); % altura
        netas(k,2)=norm(T(k+1,:)-Qo)/H2;
        if strcmp(gravitational,'yes')
            if flagbedge(k+1,1)~=0
                % para vertice no contorno do dominio
                if flagbedge(k+1,2)>200
                    gaux(k,4:6)=(R*(Qo-T(k+1,:))')';
                else
                    gaux(k,4:6)=(R*(Qo-T(k+1,:))')';
                end
            else
                % para vertice no interior do dominio
                gaux(k,4:6)=(R*(Qo-T(k+1,:))')';
            end
        end
    end

    %%%%%%%Fim do Preenchimento da segunda coluna do vetor "netas".%%%%%%%


    %Preenchimento da primeira coluna do vetor "netas".%

    v1=O(k,:)-Qo;
    v2=P(k,:)-Qo;
    ce=cross(v1,v2); % produto vetorial
    H1=norm(ce)/norm(v2); % altura
    netas(k,1)=norm(T(k,:)-Qo)/H1;
    if strcmp(gravitational,'yes')
        if flagbedge(k,1)~=0
            % para vertice no contorno do dominio
            if  flagbedge(k,2)>200

                gaux(k,1:3)=(R*(T(k,:)-Qo)')';
            else
                gaux(k,1:3)=(R*(T(k,:)-Qo)')';
            end
        else
            % para vertice no interior do dominio
            gaux(k,1:3)=(R*(T(k,:)-Qo)')';
        end
    end
    %%%%%Fim do Preenchimento da primeira coluna do vetor "netas".%%%%%%%

end

