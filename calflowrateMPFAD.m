%funçao que calcula os fluxos nas arestas internas
%equacoes 28 e 29 (heterogeneo) ou 15 e 16 (homogeneo)

function [flowrate, flowresult]=calflowrateMPFAD(p,w,s,Kde,Ded,Kn,Kt,Hesq,...
    nflagno,mobility,gravresult,gravrate,pinterp,gravno,gravelem)

global coord esurn1 esurn2 bedge inedge centelem bcflag gravitational strategy


%Initialize "bedgesize" and "inedgesize"
bedgesize = size(bedge,1);
inedgesize = size(inedge,1);

%Initialize "flowrate" and "flowresult"
flowrate = zeros(bedgesize + inedgesize,1);
flowresult = zeros(size(centelem,1),1);
m=0;
for ifacont=1:size(bedge,1);
    lef=bedge(ifacont,3);
    O=centelem(lef,:); % baricentro do elemento a esuqerda
    B1=bedge(ifacont,1);
    B2=bedge(ifacont,2);
    nor=norm(coord(B1,:)-coord(B2,:));
    
    v0=coord(bedge(ifacont,2),:)-coord(bedge(ifacont,1),:); %fase.
    v1=centelem(bedge(ifacont,3),:)-coord(bedge(ifacont,1),:);
    v2=centelem(bedge(ifacont,3),:)-coord(bedge(ifacont,2),:);
    
    if bedge(ifacont,5)<200 % se os nós esteverem na fronteira de DIRICHLET
        
        c1=nflagno(B1,2);
        c2=nflagno(B2,2);
        
        A=(Kn(ifacont)/(Hesq(ifacont)*nor));
        if strcmp(gravitational,'yes')
            if strcmp(strategy,'starnoni')
                m=gravrate(ifacont);
            else
            g1=gravno(B1,1); % gravidade no vertice 1
            g2=gravno(B2,1); % gravidade no vertice 2
            m=-A*(dot(v2,-v0)*g1+dot(v1,v0)*g2-norm(v0)^2*gravelem(lef))-...
                (g2-g1)*Kt(ifacont);
            end
        end
        auxflowrate=-A*(((O-coord(B2,:)))*(coord(B1,:)-coord(B2,:))'*c1+...
            (O-coord(B1,:))*(coord(B2,:)-coord(B1,:))'*c2-(nor^2)*p(lef))...
            -(c2-c1)*Kt(ifacont);
		
        flowrate(ifacont)=auxflowrate+m;
    else
        x=bcflag(:,1)==bedge(ifacont,5);
        r=find(x==1);
        flowrate(ifacont)= nor*bcflag(r,2);
    end
    %Attribute the flow rate to "flowresult"
    %On the left:
    flowresult(lef) = flowresult(lef) + flowrate(ifacont);
end

for iface=1:size(inedge,1)
    lef=inedge(iface,3); %indice do elemento a direita da aresta i
    rel=inedge(iface,4); %indice do elemento a esquerda da aresta i
   
    p1=pinterp(inedge(iface,1),1);
    p2=pinterp(inedge(iface,2),1);
    %calculo das vazões
        
    if strcmp(gravitational,'yes')
        if strcmp(strategy,'starnoni')
            m=gravrate(size(bedge,1)+iface,1);
        elseif strcmp(strategy,'inhouse')
            no1=inedge(iface,1);
            no2=inedge(iface,2);
            nec1=esurn2(no1+1)-esurn2(no1);
            nec2=esurn2(no2+1)-esurn2(no2);
            g1=0;
            if nflagno(no1,1)<200
                g1=gravno(no1,1);
            else
                for j=1:nec1
                    element1=esurn1(esurn2(no1)+j);
                    g1=g1+w(esurn2(no1)+j)*gravelem(element1);
                end
            end
            g2=0;
            if nflagno(no2,1)<200
                g2=gravno(no2,1);
            else
                for jj=1:nec2
                    element2=esurn1(esurn2(no2)+jj);
                    g2=g2+w(esurn2(no2)+jj)*gravelem(element2);
                end
            end
            m= Kde(iface)*(gravelem(rel,1)-gravelem(lef,1)-Ded(iface)*(g2-g1));
        end
    end
    
    
    flowrate(iface+size(bedge,1))=Kde(iface)*(p(rel)-p(lef)-Ded(iface)*(p2-p1))+m;
    %Attribute the flow rate to "flowresult"
    %On the left:
    flowresult(lef) = flowresult(lef) + flowrate(bedgesize + iface);  
    %On the right:
    flowresult(rel) = flowresult(rel) - flowrate(bedgesize + iface);
end

end