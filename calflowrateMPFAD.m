%fun�ao que calcula os fluxos nas arestas internas
%equacoes 28 e 29 (heterogeneo) ou 15 e 16 (homogeneo)

function [flowrate, flowresult]=calflowrateMPFAD(p,w,s,Kde,Ded,Kn,Kt,Hesq,...
    nflagno,mobility,gravresult,gravrate,pinterp,gravno,gravelem,grav_elem_escalar)

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

    if bedge(ifacont,5)<200 % se os n�s esteverem na fronteira de DIRICHLET

        c1=nflagno(B1,2);
        c2=nflagno(B2,2);

        A=(Kn(ifacont)/(Hesq(ifacont)*nor));
        if strcmp(gravitational,'yes')
            m=gravrate(ifacont);

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
    %calculo das vaz�es

    if strcmp(gravitational,'yes')

        m=gravrate(size(bedge,1)+iface,1);
    end
    flowrate(iface+size(bedge,1))=Kde(iface)*(p(rel)-p(lef)-Ded(iface)*(p2-p1))+m;
    %Attribute the flow rate to "flowresult"
    %On the left:
    flowresult(lef) = flowresult(lef) + flowrate(bedgesize + iface);
    %On the right:
    flowresult(rel) = flowresult(rel) - flowrate(bedgesize + iface);
end

end