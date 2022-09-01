
function [x, i,ciclos,tolerancia]=broyden(M_old,RHS_old,p_old,tol,kmap,parameter,...
   w,s,nflagface,fonte,gamma,nflagno,p_old1,...
    weightDMP,auxface,calnormface,wells,mobility,gravresult,varargin)


% Parse any name-value pairs that will override the defaults:
p = inputParser;
%substituir por tolnewton
addParameter(p,'tolfun',1e-5);  % Default values
%substituir por tolnewton
addParameter(p,'tolx',1e-4);
%valor qualquer pode ser escolhido (broyden original não limita por
%iteração)
addParameter(p,'maxiter',30);
addParameter(p,'bounds',[]);
%

parse(p,varargin{:})

% Fields of opt take the default values unless overriden by a
% name/value pair


x = p_old;  %%P_OLD
f = M_old*x-RHS_old;
R0 = norm(f); %%CÁLCULO DE F(X0) PARA UTILIZAÇÃO E MEDIÇÃO DE UM ERRO APRESENTADO
if ~(size(x) == size(f))
    error('f must return a column vector of the same size as x0')
end
%calculo da matriz jacobiana
%J = jacobi(fun,x);  % Intial Jacobian matrix %%M_OLD E P_OLD

J=aproxmjacobian(f,p_old1,p_old,nflagface,nflagno,w,s,...
    parameter,weightDMP,kmap,fonte,auxface,mobility,0, 0, 0, 0, ...
    0,calnormface,wells);

ciclos=1;
i=1;
er=1;
while  tol<er
    if any(isnan(J(:))) || rcond(full(J)) < 1e-15
        error('Singular jacobian at iteration %d\n',i)
    end
    
    dx=-J\f; %valor da jacobiana 
    
    x  = x+dx; % perturbação
    % calculo da interpolação
    [pinterp_new]=pressureinterp(x,nflagface,nflagno,w,s,...
        parameter,weightDMP);
    % calculo da matriz 
    [M_new,RHS_new]=globalmatrix(x,pinterp_new,gamma,nflagface,nflagno...
        ,parameter,kmap,fonte,metodoP,w,s,weightDMP,auxface,...
        wells,mobility,0, 0, 0, 0, 0,calnormface);
    f=M_new*x-RHS_new;
    
    J  = J + f*dx'/(dx'*dx);
    R = norm(f);
    
    if (R0 ~= 0.0)
        er = abs(R/R0)
    else
        er = 0;
    end
    
    i=i+1;
end
tolerancia= er;
end

