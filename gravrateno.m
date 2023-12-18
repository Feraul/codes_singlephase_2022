function [g]=gravrateno(No,gravrate,V)
global bedge inedge esurn2 esurn1 

% faces na vizinhança do No
facesV=V(1:2,:,No);
m=1;
for j=1:size(facesV,2)
    % elemento atual
    elemaux=esurn1(esurn2(No)+j);
    if max(facesV(:,j))~=0
        for i=1:size(facesV,1)
            if facesV(i,j)<size(bedge,1) || facesV(i,j)==size(bedge,1)
               % aloca os fluxos nas interfaces
                g(m,i)=gravrate(facesV(i,j));
            else
               ifaceaux= facesV(i,j)-size(bedge,1);
               lef= inedge(ifaceaux,3);
               if elemaux==lef
                   % fluxo sai com sinal positivo
                  g(m,i)=gravrate(facesV(i,j));
               else
                  g(m,i)=-gravrate(facesV(i,j));
               end
            end
        end
        m=m+1;
    end
end

end