function [gravaux]=gravinterpHP(weightDMP,parameter,gravface,grav_elem_escalar,nflagface)
global inedge coord bedge bcflag 

    mobil=1;
    mobilO=1;
    gravaux=zeros(size(inedge,1)+size(bedge,1),1);
    %% interpola��o das press�es no pontos armonicos internos
    gravaux(size(bedge,1)+1:size(bedge,1)+size(inedge,1),1)=weightDMP(:,1).*grav_elem_escalar(weightDMP(:,3))+ weightDMP(:,2).*grav_elem_escalar(weightDMP(:,4));
    %% interpola��o das press�es no pontos medios do contorno
    for ifacont=1:size(bedge,1)
        lef=bedge(ifacont,3);
        
        % Quando o ponto harmonico oposto pertece a face interior da malha ou
        % contorno de Dirichlet
        
        if bedge(ifacont,5)>200
            auxfacelef1=parameter(1,3,ifacont);
            auxfacelef2=parameter(1,4,ifacont);
            if auxfacelef1==ifacont
                faceoposto=auxfacelef2;
                atualksi=parameter(1,1,ifacont);
                opostoksi=parameter(1,2,ifacont);
            else
                faceoposto=auxfacelef1;
                atualksi=parameter(1,2,ifacont);
                opostoksi=parameter(1,1,ifacont);
            end
            normcont=norm(coord(bedge(ifacont,1),:)-coord(bedge(ifacont,2),:));
            x=bcflag(:,1)==bedge(ifacont,5);
            r=find(x==1);
            
            % calcula o fluxo na face "ifacelef2"
            fluxoN=normcont*bcflag(r,2);
            %mobil=mobility(ifacont);
            % "faceoposto" � aquele face oposto ao "ifacont"
            if faceoposto<size(bedge,1) || faceoposto==size(bedge,1)
                % faceoposto pode ser contorno de Dirichlet
                if nflagface(faceoposto,1)<200
                    %  caso I. Quando a face oposto pertence ao contorno de Dirichlet
                    pressure=gravface(faceoposto);
                    
                    termo1=(opostoksi+atualksi)/atualksi;
                    termo2= fluxoN/(mobil*normcont*atualksi);
                    termo3= opostoksi/atualksi;
                    
                    gravaux(ifacont,1)= termo1*grav_elem_escalar(lef)-termo2-termo3*pressure;
                else
                    %  caso II. Quando a face oposto pertence ao contorno de Neumann
                    %mobilO=mobility(faceoposto);
                    normcontopost=norm(coord(bedge(faceoposto,1),:)-coord(bedge(faceoposto,2),:));
                    
                    x=bcflag(:,1)==bedge(faceoposto,5);
                    r=find(x==1);
                    fluxOpost=normcontopost*bcflag(r,2);
                    
                    if auxfacelef1==faceoposto
                        atualksO=parameter(1,1,faceoposto);
                        opostksO=parameter(1,2,faceoposto);
                    else
                        atualksO=parameter(1,2,faceoposto);
                        opostksO=parameter(1,1,faceoposto);
                    end
                    
                    sumaksiatual=  opostksO*(atualksi + opostoksi);
                    
                    sumaksiopost= opostoksi*(atualksO + opostksO);
                    
                    
                    sumatotalnumer=  sumaksiatual - sumaksiopost;
                    
                    sumatotaldenom= atualksi*opostksO - atualksO*opostoksi;
                    
                    termo1=sumatotalnumer/sumatotaldenom;
                    
                    if abs(sumatotaldenom)<1e-5
                        termo2=0;
                    else
                        termo2=(opostoksi*(fluxOpost/(mobilO*normcontopost)) - opostksO*(fluxoN/(mobil*normcont)))/sumatotaldenom;
                    end
                    gravaux(ifacont,1)= termo1*grav_elem_escalar(lef)+termo2;
                    
                end
            else
                % caso III. Quando a face oposto pertence ao interior da
                % malha
                termo1=(opostoksi+atualksi)/atualksi;
                
                termo2=(fluxoN/(mobil*normcont*atualksi));
                
                % Calculo das contribui��es do elemento a esquerda
                auxweightlef1=weightDMP(faceoposto-size(bedge,1),1);
                auxweightlef2=weightDMP(faceoposto-size(bedge,1),2);
                
                auxlef=weightDMP(faceoposto-size(bedge,1),3);
                auxrel=weightDMP(faceoposto-size(bedge,1),4);
                if auxlef==lef
                    auxelematual=auxlef;
                    auxelemopost=auxrel;
                    pesatual=auxweightlef1;
                    pesopost=auxweightlef2;
                else
                    auxelematual=auxrel;
                    pesatual=auxweightlef2;
                    pesopost=auxweightlef1;
                    auxelemopost=auxlef;
                end
                termo3= (opostoksi/atualksi);
                
                pfaceopost= pesatual*grav_elem_escalar(auxelematual)+pesopost*grav_elem_escalar(auxelemopost);
                
                gravaux(ifacont,1)= termo1*grav_elem_escalar(auxelematual)-termo2-termo3*pfaceopost;
            end
        else
            gravaux(ifacont,1)= gravface(ifacont); 
        end
    end
end