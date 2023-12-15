% objetivo: Montagem da matriz global M e I
function [ M, I ] = globalmatrixmpfad( w,s, Kde, Ded, Kn, Kt, nflagno, ...
    Hesq,fonte,gravresult,gravrate,N)

global coord elem esurn1 esurn2  bedge inedge  centelem bcflag gravitational...
    strategy elemarea

%-----------------------inicio da rutina ----------------------------------%
%Constroi a matriz global.

M=sparse(size(elem,1),size(elem,1)); %Prealocacao de M.
I=sparse(size(elem,1),1);
% contribuicao do termo de fonte
I=I+fonte;
grav_no1=0;
grav_no2=0;
m=0;

for ifacont=1:size(bedge,1)
    lef=bedge(ifacont,3);
    
    v0=coord(bedge(ifacont,2),:)-coord(bedge(ifacont,1),:); %fase.
    v1=centelem(bedge(ifacont,3),:)-coord(bedge(ifacont,1),:);
    v2=centelem(bedge(ifacont,3),:)-coord(bedge(ifacont,2),:);
    normcont=norm(v0);
    
    % Tratamento do nó nos vértices 2 e 4%
    A=-Kn(ifacont)/(Hesq(ifacont)*norm(v0));
    
    if bedge(ifacont,5)<200
        % Contorno de Dirichlet
        c1=nflagno(bedge(ifacont,1),2);
        c2=nflagno(bedge(ifacont,2),2);
        
        %Preenchimento do termo gravitacional
        
        if strcmp(gravitational,'yes')
            if strcmp(strategy,'starnoni')||strcmp(strategy,'inhouse')|| strcmp(strategy,'GravConsist')
                m=gravrate(ifacont,1);
            end
        end
        %montagem da matriz global
        M(lef,lef)=M(lef,lef)-A*(norm(v0)^2);
        % termo de fonte
        I(lef)=I(lef)-A*(dot(v2,-v0)*c1+dot(v1,v0)*c2)+(c2-c1)*Kt(ifacont)+m;
    else
        
        if strcmp(gravitational,'yes')
            if strcmp(strategy,'starnoni')||strcmp(strategy,'inhouse')|| strcmp(strategy,'GravConsist')
                m=gravrate(ifacont,1);
            end
        end
        % Contorno de Neumann
        x=bcflag(:,1)==bedge(ifacont,5);
        r=find(x==1);
        I(lef)=I(lef) -normcont*bcflag(r,2)+m;
    end
end

% contribuição nas faces internas
for iface=1:size(inedge,1)
    lef=inedge(iface,3);
    rel=inedge(iface,4);
    
    %Contabiliza as contribuições do fluxo numa aresta para os elementos %
    %a direita e a esquerda dela.                                        %
    
    M(lef, lef)=M(lef, lef)- Kde(iface);
    M(lef, rel)=M(lef, rel)+ Kde(iface);
    M(rel, rel)=M(rel, rel)- Kde(iface);
    M(rel, lef)=M(rel, lef)+ Kde(iface);
    
    %Se os nós das arestas estiverem em fronteiras de Dirichlet, suas
    %contribuições serão contabilizadas logo abaixo.
    
    if nflagno(inedge(iface,1),1)<200
        I(lef)=I(lef)-Kde(iface)*Ded(iface)*nflagno(inedge(iface,1),2);
        I(rel)=I(rel)+Kde(iface)*Ded(iface)*nflagno(inedge(iface,1),2);
    end
    if nflagno(inedge(iface,2),1)<200
        I(lef)=I(lef)+Kde(iface)*Ded(iface)*nflagno(inedge(iface,2),2);
        I(rel)=I(rel)-Kde(iface)*Ded(iface)*nflagno(inedge(iface,2),2);
    end
    % quando o nó pertece ao contorno de Neumann
    if nflagno(inedge(iface,1),1)==201
        
        I(inedge(iface,3))=I(inedge(iface,3))-Kde(iface)*Ded(iface)*s(inedge(iface,1)); %ok
        I(inedge(iface,4))=I(inedge(iface,4))+Kde(iface)*Ded(iface)*s(inedge(iface,1)); %ok
    end
    if nflagno(inedge(iface,2),1)==201
        
        I(inedge(iface,3))=I(inedge(iface,3))+Kde(iface)*Ded(iface)*s(inedge(iface,2)); %ok
        I(inedge(iface,4))=I(inedge(iface,4))-Kde(iface)*Ded(iface)*s(inedge(iface,2)); %ok
        
    end
    
    %Contabilização das contribuições dos nós que não estão na
    %fronteiras de Dirichlet.
    % first node
    if nflagno(inedge(iface,1),1)>200
        for j=1:(esurn2(inedge(iface,1)+1)-esurn2(inedge(iface,1)))
            
            post_cont=esurn2(inedge(iface,1))+j;
            
            M(lef, esurn1(post_cont))=M(lef,esurn1(post_cont)) + Kde(iface)*Ded(iface)*w(post_cont);
            M(rel, esurn1(post_cont))=M(rel,esurn1(post_cont)) - Kde(iface)*Ded(iface)*w(post_cont);
            
        end
        
    end
    % second node
    if nflagno(inedge(iface,2),1)>200
        for j=1:(esurn2(inedge(iface,2)+1)-esurn2(inedge(iface,2)))
            post_cont=esurn2(inedge(iface,2))+j;
            
            M(lef, esurn1(post_cont))=M(lef,esurn1(post_cont)) - Kde(iface)*Ded(iface)*w(post_cont);
            M(rel, esurn1(post_cont))=M(rel,esurn1(post_cont)) + Kde(iface)*Ded(iface)*w(post_cont);
        end
    end
    
    % termo gravitacional
    if strcmp(gravitational,'yes')
        m=gravrate(size(bedge,1)+iface,1);
        I(lef)=I(lef)+m ;
        I(rel)=I(rel)-m ;
    end
end
end

