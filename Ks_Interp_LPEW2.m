function [ Kt1, Kt2, Kn1, Kn2,gaux3 ] = Ks_Interp_LPEW2( O, T, Qo, kmap, No,gravelem)

global esurn2 esurn1 elem  gravitational
%Retorna os K(n ou t) necess�rios para a obten��o dos weights. kmap � a
%matriz de permeabilidade; Ex: Kt1->linhaN=Kt1(cellN);

nec=esurn2(No+1)-esurn2(No);

Kt1=zeros(nec,2); %As colunas representam i=1 e i=2.
Kt2=zeros(nec,1);
Kn1=zeros(nec,2);
Kn2=zeros(nec,1);
K=zeros(3);
R=[0 1 0; -1 0 0; 0 0 0];
gaux3=0;
%Constru��o do tensor permeabilidade.%

%C�lculo das primeiras constantes, para todas as c�lulas que concorrem num
%vertice "No".
for k=1:nec
    % elemento j
    j=esurn1(esurn2(No)+k);
    % permeabilidade
    K(1,1)=kmap(elem(j,5),2);
    K(1,2)=kmap(elem(j,5),3);
    K(2,1)=kmap(elem(j,5),4);
    K(2,2)=kmap(elem(j,5),5);
    for i=1:2
        if (size(T,1)==size(O,1))&&(k==nec)&&(i==2)
            Kn1(k,i)=((R*(T(1,:)-Qo)')'*K*(R*(T(1,:)-Qo)'))/norm(T(1,:)-Qo)^2;
            Kt1(k,i)=((R*(T(1,:)-Qo)')'*K*(T(1,:)-Qo)')/norm(T(1,:)-Qo)^2;
        else
            Kn1(k,i)=((R*(T(k+i-1,:)-Qo)')'*K*(R*(T(k+i-1,:)-Qo)'))/norm(T(k+i-1,:)-Qo)^2;
            Kt1(k,i)=((R*(T(k+i-1,:)-Qo)')'*K*(T(k+i-1,:)-Qo)')/norm(T(k+i-1,:)-Qo)^2;
        end
    end
    %------------------------- Tensores ----------------------------------%
    if (size(T,1)==size(O,1))&&(k==nec)
        %------------ Calculo dos K's internos no elemento ---------------%
        Kn2(k)=((R*(T(1,:)-T(k,:))')'*K*(R*(T(1,:)-T(k,:))'))/norm(T(1,:)-T(k,:))^2;
        Kt2(k)=((R*(T(1,:)-T(k,:))')'*K*(T(1,:)-T(k,:))')/norm(T(1,:)-T(k,:))^2;
        vec= (R*(T(1,:)-T(k,:))')';
    else
        Kn2(k)=(R*(T(k+1,:)-T(k,:))')'*K*(R*(T(k+1,:)-T(k,:))')/norm(T(k+1,:)-T(k,:))^2;
        Kt2(k)=((R*(T(k+1,:)-T(k,:))')'*K*(T(k+1,:)-T(k,:))')/norm(T(k+1,:)-T(k,:))^2;
        vec=(R*(T(k+1,:)-T(k,:))')';
    end
    if strcmp(gravitational,'yes')
        gaux3(k)=dot(vec*K,gravelem(j,:));
    end
end
end

