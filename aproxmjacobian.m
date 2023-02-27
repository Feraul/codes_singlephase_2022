function [J]=aproxmjacobian(Fk,p_new,p_old,nflagface,nflagno,w,s,...
    parameter,weightDMP,kmap,fonte,auxface,mobility,Hesq, Kde, Kn, Kt, ...
    Ded,calnormface,wells,M_old)

global elem
nelem=size(elem,1);
J=sparse(nelem,nelem);
x=p_old;
I=eye(nelem);
MMM=full(M_old)~=0;
delta = 1e-3*sqrt(norm(p_old));

for ielem=1:nelem
    %  xi+h
    x=x+ delta*MMM(:,ielem);
%     x=[1.0026
%     1.0026
%     1.0000
%     1.0026
%     1.0000
%     1.0036
%     1.0036
%     1.0036
%     1.0000
%     1.0036
%     1.0036
%     1.0036
%     1.0000
%     1.0036
%     1.0000
%     1.0036];
    % interpolation point
    [pinterp_new]=pressureinterp(x,nflagface,nflagno,w,s,...
        parameter,weightDMP,mobility);
    
    % Calculo da matriz global
    [auxM,auxRHS]=globalmatrix(x,pinterp_new,0,nflagface,nflagno...
        ,parameter,kmap,fonte,w,s,weightDMP,wells,...
        mobility,Hesq, Kde, Kn, Kt, Ded,calnormface,0,0,0,0,0,0);
    
    % f(xk+1)
    Fkk= auxM*x - auxRHS;
    
    % Montagem da matriz Jacobiano por método Diferencia Finita
    J(1:nelem,ielem)=(Fkk(:)-Fk(:))./(delta);
    
    % Atualiza o vetor "x"
    x=p_old;
    delta=1e-3*sqrt(norm(x));
end

end