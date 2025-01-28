function [ w,s,wg] = Pre_LPEW_2(kmap,N,gravrate,gradgravelem,V)
global coord bcflag bedge nsurn1 nsurn2 gravitational esurn1 esurn2 elem
% devolve os pesos "w" cujos elementos são organizados em um vetor
%Retorna todos os parâmetros necessários às expressões dos fluxos.

apw=ones(size(coord,1),1);
apw1=ones(size(coord,1),1);
r=zeros(1,2);
s=zeros(size(bedge,1),1);
K=zeros(3);
for No=1:size(coord,1)
    
    %[g]=gravnode(N,kmap,No,gravelem);
    %     if strcmp(gravitational,'yes')
    %         [gaux1]=gravrateno(No,gravrate,V);
    %     end
    % calcula
    % O--> coordenadas do baricentro na vizinhança do nó "No"
    % P--> coordenadas dos vértices na vizinhança do nó "No"
    % T--> coordenadas dos pontos medios nas vizinhas ao nó "No"
    % Qo-> coordenada do nó em questão
    [ O, P, T, Qo,flagbedge ] = OPT_Interp_LPEW(No);
    % calcula os angulos apropiados para calculas os pesos
    
    [ ve2, ve1, theta2, theta1 ] = angulos_Interp_LPEW2( O, P, T, Qo,No );
    % calculas as netas uma relacao de de alturas
    
    [ neta,gaux,gauxm] = netas_Interp_LPEW( O, P, T, Qo, No,kmap,gradgravelem,flagbedge);
    % calculas as projecoes normais em torno do nó "No"
    [ Kt1, Kt2, Kn1, Kn2,gaux3 ] = Ks_Interp_LPEW2( O, T, Qo, kmap, No,gradgravelem);
    % calcula os lambdas
    [ lambda,r,gaux2 ] =  Lamdas_Weights_LPEW2( Kt1, Kt2, Kn1, Kn2, theta1,...
        theta2, ve1, ve2, neta, P, O,r,gaux);
   % [ lambda,r,gaux2m ] =  aux_Lamdas_Weights_LPEW2( Kt1, Kt2, Kn1, Kn2, theta1,...
   %    theta2, ve1, ve2, neta, P, O,r,gauxm);
   
    for k=0:size(O,1)-1
        w(apw(No)+k,1)=lambda(k+1)/sum(lambda); %calculo dos pesos 
      
    end
    
    apw(No+1)=apw(No)+size(O,1);
    c=esurn2(No+1)-esurn2(No);
    m=0;
    for i=1:c
        j=esurn1(esurn2(No)+i);
        %Essa é UMA maneira de construir os tensores
        K(1,1) = kmap(elem(j,5),2);
        K(1,2) = kmap(elem(j,5),3);
        K(2,1) = kmap(elem(j,5),4);
        K(2,2) = kmap(elem(j,5),5);
        
        %gaux2: equacao 4.33 do manual 
        %m(i)=dot((K*gradgravelem(j,:)')',gaux2(i,:));
        m=m+dot((K*gradgravelem(j,:)')',gaux2(i,:));
    end
    %for k=0:size(O,1)-1
    %    wg(apw(No)+k,1)=(m(k+1))/sum(lambda); %calculo dos pesos  
    %end
    %wg(No,1)=m/sum(lambda); % erro -> 0.0023
    wg(No,1)=(sum(gaux3)+m)/sum(lambda); % erro -> 0.0038
    %wg(No,1)=sum(gaux2m)/sum(lambda); % erro-> 0.0035
    %wg(No,1)=(sum(gaux3)+sum(gaux2m))/sum(lambda); % erro -> 0.005
    %wg(No,1)=(sum(ww.*gaux3))/sum(lambda); % erro -> 0.0054
        
    % calculando os pesos nos vertices do contorno de Neumann
    % N ordena faces na vizinhanca de um vertices, comecando pela face do
    % contorno
    vetor = nsurn1(nsurn2(No) + 1:nsurn2(No + 1));
    comp1 = N(No,1);
    comp2 = N(No,length(vetor));
    MM=bedge(:,1)==No;
    MMM= find(MM == 1);
    % Se o vertices 'No' pertence ao contorno de Neumann
    if comp1<= size(bedge,1) && comp2 <=size(bedge,1) && 200<bedge(MMM,4)
        % norma das faces
        norm1=norm(coord(bedge(comp1,2),1:2)-coord(bedge(comp1,1),1:2));
        norm2=norm(coord(bedge(comp2,2),1:2)-coord(bedge(comp2,1),1:2));
        a = bcflag(:,1) == bedge(comp1,5);
        s1 = find(a == 1);
        b = bcflag(:,1) == bedge(comp2,5);
        s2 = find(b == 1);
        %------------------------------------------------------------------
        if strcmp(gravitational,'yes')
             % contribuicoes do termo gravitacional
           R=[0 1; -1 0 ];
            %------------------------------
            K1(1,1)=kmap(elem(bedge(comp1,3),5),2);
            K1(1,2)=kmap(elem(bedge(comp1,3),5),3);
            K1(2,1)=kmap(elem(bedge(comp1,3),5),4);
            K1(2,2)=kmap(elem(bedge(comp1,3),5),5);
            
            v1=0.5*(coord(bedge(comp1,2),1:2)+coord(bedge(comp1,1),1:2))-coord(No,1:2);
            %------------------------------
            K2(1,1)=kmap(elem(bedge(comp2,3),5),2);
            K2(1,2)=kmap(elem(bedge(comp2,3),5),3);
            K2(2,1)=kmap(elem(bedge(comp2,3),5),4);
            K2(2,2)=kmap(elem(bedge(comp2,3),5),5);

            v2=0.5*(coord(bedge(comp2,2),1:2)+coord(bedge(comp2,1),1:2))-coord(No,1:2);

            %-------------------------------
            m1=dot((R*(v1)')'*K1,gradgravelem(bedge(comp1,3),1:2));
            m2= dot((R*(v2)')'*K2,gradgravelem(bedge(comp2,3),1:2));
           
        else
            m1=0;
            m2=0;
        end
        % da errado quando colocamos o termo gravitacional
        %s(No,1) = -(1/sum(lambda))*(r(1,1)*(norm1*bcflag(s1,2))+...
        %                            r(1,2)*(norm2*bcflag(s2,2)));

        %wg(No,1)=wg(No,1)+(1/sum(lambda))*(r(1,1)*m1+r(1,2)*m2); 
        %esta rutina funciona quando o vertice da quina da malha
        %computacional pertence ao contorno de Dirichlet
        s(No,1) = -(1/sum(lambda))*(r(1,1)*bcflag(s1,2) + ...
            r(1,2)*bcflag(s2,2));
        %esta rutina funciona quando o vertice da quina da malha
        %computacional pertence ao contorno de Dirichlet ou Neumann
        %s(No,1) = -(1/sum(lambda))*(r(No,1)*m11 + ...
        %    r(No,2)*m12);
    end  %End of IF
end
end

