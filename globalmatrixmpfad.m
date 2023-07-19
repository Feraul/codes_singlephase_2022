% objetivo: Montagem da matriz global M e I
function [ M, I ] = globalmatrixmpfad( w,s, Kde, Ded, Kn, Kt, nflagno, ...
    Hesq,fonte,gravresult,gravrate,gravno,gravelem,grav_elem_escalar,wg)

global coord elem esurn1 esurn2  bedge inedge  centelem bcflag gravitational...
    strategy elemarea

%-----------------------inicio da rutina ----------------------------------%
%Constroi a matriz global.

M=sparse(size(elem,1),size(elem,1)); %Prealocacao de M.
I=sparse(size(elem,1),1);
% contribuicao do termo de fonte
I=I+fonte;

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
            if strcmp(strategy,'starnoni')||strcmp(strategy,'inhouse3')
                m=gravrate(ifacont);
            elseif strcmp(strategy,'inhouse1')
                g1=gravno(bedge(ifacont,1),1); % gravidade no vertice 1
                g2=gravno(bedge(ifacont,2),1); % gravidade no vertice 2
               % m=-(A*(dot(v2,-v0)*g1+dot(v1,v0)*g2-(norm(v0)^2*grav_elem_escalar(lef)))-(g2-g1)*Kt(ifacont));
                m=gravrate(ifacont)+g2+g1;
            end
        else
            m=0;
        end
        % ambos os vertices pertenecem ao contorno de Dirichlet
        if nflagno(bedge(ifacont,2),1)<200 && nflagno(bedge(ifacont,1),1)<200
            %montagem da matriz global 
            M(lef,lef)=M(lef,lef)-A*(norm(v0)^2);
            % termo de fonte
            I(lef)=I(lef)-A*(dot(v2,-v0)*c1+dot(v1,v0)*c2)+(c2-c1)*Kt(ifacont)+m;
        else
            % quando um dos vertices da quina da malha computacional
            % pertence ao contorno de Neumann
            if nflagno(bedge(ifacont,1),1)>200
                %montagem da matriz global
                M(lef,lef)=M(lef,lef)-A*(norm(v0)^2)+Kt(ifacont)+A*dot(v2,-v0);
                % termo de fonte
                I(lef)=I(lef)-A*(dot(v1,v0)*c2)+(c2)*Kt(ifacont)+m;
            elseif nflagno(bedge(ifacont,2),1)>200
                %montagem da matriz global
                M(lef,lef)=M(lef,lef)-A*(norm(v0)^2)-Kt(ifacont)+A*dot(v1,v0);
                % termo de fonte
                I(lef)=I(lef)-A*(dot(v2,-v0)*c1)+(-c1)*Kt(ifacont)+m;
            end
        end
    else
        % Contorno de Neumann
        x=bcflag(:,1)==bedge(ifacont,5);
        r=find(x==1);
        I(lef)=I(lef) -normcont*bcflag(r,2);
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
        %I(lef)=I(lef)+  wg(inedge(iface,1));
        %I(rel)=I(rel)-   wg(inedge(iface,1));
    end
    % second node
    if nflagno(inedge(iface,2),1)>200
        for j=1:(esurn2(inedge(iface,2)+1)-esurn2(inedge(iface,2)))
            
            post_cont=esurn2(inedge(iface,2))+j;
            
            M(lef, esurn1(post_cont))=M(lef,esurn1(post_cont)) - Kde(iface)*Ded(iface)*w(post_cont);
            
            M(rel, esurn1(post_cont))=M(rel,esurn1(post_cont)) + Kde(iface)*Ded(iface)*w(post_cont);
            
        end
        %I(lef)=I(lef)-  wg(inedge(iface,2));
        %I(rel)=I(rel)+  wg(inedge(iface,2));
    end
    % termo gravitacional
    if strcmp(gravitational,'yes')
        if strcmp(strategy,'starnoni') ||strcmp(strategy,'inhouse3')
            m=gravrate(size(bedge,1)+iface,1);
        elseif strcmp(strategy,'inhouse')
            no1=inedge(iface,1);
            no2=inedge(iface,2);
            %g1=gravno(no1,1);
            %g2=gravno(no2,1);
            g1=0;
            nec1=esurn2(no1+1)- esurn2(no1);
            nec2=esurn2(no2+1)- esurn2(no2);
            for j=1:nec1
                    element1=esurn1(esurn2(no1)+j);
                    g1=g1+w(esurn2(no1)+j)*grav_elem_escalar(element1);
            end
            g2=0;
            for j=1:nec2
                    element2=esurn1(esurn2(no2)+j);
                    g2=g2+w(esurn2(no2)+j)*grav_elem_escalar(element2);
            end
           % m= -Kde(iface)*(grav_elem_escalar(rel)-grav_elem_escalar(lef)-Ded(iface)*(g2-g1));
            m= gravrate(size(bedge,1)+iface,1)+g2+g1;
        else
            m=0;
        end
        
        I(lef)=I(lef)+m ;
        I(rel)=I(rel)-m ;
    end
end
end

