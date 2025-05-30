function [p_new,tabletol,iter,ciclos]=picardMPE(M_old,RHS_old,nitpicard,tolpicard,kmap,...
    parameter,w,s,nflagface,fonte,p_old,gamma,nflagno,...
    weightDMP,auxface,wells,mobility,Hesq, Kde, Kn, Kt, Ded,calnormface)
%% calculo do residuo Inicial
R0=norm(M_old*p_old-RHS_old);
%% inicializando dados para itera��o Picard
ciclos=9;
%=======================================================================
       %% inicializando dados para itera��o Picard
% if rcond(full(M_old))<1e-5
     [L,U] = ilu(M_old,struct('type','ilutp','droptol',1e-6));
     
     [p_old,]=bicgstab(M_old,RHS_old,1e-11,1000,L,U,p_old);
%  else
%     [p_old,]=bicgstab(M_old,RHS_old,1e-11,1000);
%  end
  %[p_old,fl1,rr1,it1,rv1]=gmres(M_old,RHS_old,10,1e-11,1000,L,U,p_old);
%==========================================================================
 [p_new,iter,tabletol] = extrapolate(p_old,ciclos, 10000, 'MPE',nflagface,nflagno,w,s,...
        parameter,weightDMP,kmap,fonte,auxface,mobility,Hesq, Kde, Kn, Kt,...
        Ded,calnormface,wells,tolpicard,R0);
%=======================================================================
end