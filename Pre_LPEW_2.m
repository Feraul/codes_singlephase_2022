function [ w,s] = Pre_LPEW_2(kmap,N,gravrate,nflagface)
global coord bcflag bedge nsurn1 nsurn2 gravitational 
% devolve os pesos "w" cujos elementos s?o organizados em um vetor
%Retorna todos os par?metros necess?rios ?s express?es dos fluxos.
apw=ones(size(coord,1),1);
r=zeros(size(coord,1),2);
for No=1:size(coord,1)
    % calcula
    % O--> coordenadas do baricentro na vizinhan?a do n? "No"
    % P--> coordenadas dos v?rtices na vizinhan?a do n? "No"
    % T--> coordenadas dos pontos medios nas vizinhas ao n? "No"
    % Qo-> coordenada do n? em quest?o
    [ O, P, T, Qo ] = OPT_Interp_LPEW(No);
    % calcula os angulos apropiados para calculas os pesos
    
    [ ve2, ve1, theta2, theta1 ] = angulos_Interp_LPEW2( O, P, T, Qo,No );
    % calculas as netas uma relacao de de alturas
    
    [ neta ] = netas_Interp_LPEW( O, P, T, Qo, No );
    % calculas as projecoes normais em torno do n? "No"
    [ Kt1, Kt2, Kn1, Kn2 ] = Ks_Interp_LPEW2( O, T, Qo, kmap, No);
    % calcula os lambdas
    [ lambda,r ] =  Lamdas_Weights_LPEW2( Kt1, Kt2, Kn1, Kn2, theta1,...
        theta2, ve1, ve2, neta, P, O,Qo,No,T,r );
    for k=0:size(O,1)-1
        w(apw(No)+k,1)=lambda(k+1)/sum(lambda); %calculo dos pesos
    end
    
    apw(No+1)=apw(No)+size(O,1);
    % calculando os pesos na condi??o de contorno de Neumann para poder
    % interpolar as pressoes naquele contorno
    % N ordena faces na vizinhanca de um vertices, comecando pela face do
    % contorno 
    vetor = nsurn1(nsurn2(No) + 1:nsurn2(No + 1));
    comp1 = N(No,1);
    comp2 = N(No,length(vetor));
    if comp1<= size(bedge,1) && comp2 <=size(bedge,1)
        a = bcflag(:,1) == bedge(comp1,5);
        s1 = find(a == 1);
        m11=nflagface(comp1,2);
        b = bcflag(:,1) == bedge(comp2,5);
        s2 = find(b == 1);
        m12=nflagface(comp2,2);
        if strcmp(gravitational,'yes')
            m1= gravrate(comp1,1);
            m2= gravrate(comp2,1);
        else
            m1=1;
            m2=1;
        end
        % da errado quando colocamos o termo gravitacional
%        s(No,1) = -(1/sum(lambda))*(r(No,1)*m11 + r(No,1)*m1+ ...
%            r(No,2)*m12+r(No,2)*m2);
        %esta rutina funciona quando o vertice da quina da malha
        %computacional pertende ao contorno de Dirichlet
%         s(No,1) = -(1/sum(lambda))*(r(No,1)*bcflag(s1,2) + ...
%            r(No,2)*bcflag(s2,2));
        %esta rutina funciona quando o vertice da quina da malha
        %computacional pertence ao contorno de Dirichlet ou Neumann
       s(No,1) = -(1/sum(lambda))*(r(No,1)*m11 + ...
           r(No,2)*m12);
    end  %End of IF
end
end

