% este codigo contem todas as informacoes dos diferentes casos, como por
% exemplo: tensor de permeabilidade (kmap), pressão analitica (u), termos de fonte (fonte),
% velocidade analitica (vel),gravidade (grav)
function[elem,kmap,normKmap,u,bedge,fonte,vel,grav,gravno,gravface,...
    grav_elem_escalar]=benchmarks(kmap,elem,bedge)
global centelem coord inedge normals elemarea bcflag benchmark
normKmap=0;
vel=0;
u=0;
fonte=0;
grav=zeros(size(elem,1),3);
gravno=0;
gravface=0;
grav_elem_escalar=0;

switch benchmark
    case 'edwards'
        for i = 1:size(centelem,1)
            %Define "x" and "y"
            x = centelem(i,1);
            y = centelem(i,2);
            
            if x<=0.5
                u(i,1)= 10+2*x*y;
                kmap(i,:) = [i 1 0.5 0.5 1];
                fonte(i,1)=-20*elemarea(i);
            else
                u(i,1)= 10.75-1.5*x+9*y+2*x*y;
                kmap(i,:) = [i 10 2 2 100];
                fonte(i,1)=-8*elemarea(i);
            end
            elem(i,5)=i;
        end
        K=kmap;
    case 'miao'
        for i = 1:size(centelem,1)
            %Define "x" and "y"
            x = centelem(i,1);
            y = centelem(i,2);
            if x<=0.5
                u(i,1)= 14*x+y;
                kmap(i,:) = [i 3 1 1 3];
                
            else
                u(i,1)= 4*x+y+5;
                kmap(i,:) = [i 10 3 3 10];
            end
            elem(i,5)=i;
        end
        R1=[0 -1 0; 1 0 0; 0 0 0];
        for iface=1:size(bedge,1)+size(inedge,1)
            
            if iface< size(bedge,1) || iface==size(bedge,1)
                v1=bedge(iface,1);
                v2=bedge(iface,2);
            else
                v1=inedge(iface-size(bedge,1),1);
                v2=inedge(iface-size(bedge,1),2);
            end
            IJ=coord(v2,:)-coord(v1,:);
            norma=norm(IJ);
            nij=R1*IJ'/norma;
            p1=(coord(v2,:)+coord(v1,:))*0.5;
            if p1(1,1)<=0.5
                a=[-43, -17, 0];
                
            else
                a=[-43, -22, 0];
            end
            F(iface,1)=dot(a,nij');
            
        end
        vel=F;
        K=kmap;
    case 'starnonigrav1'
        
        R=[0 1 0; -1 0 0; 0 0 0];
        for i = 1:size(centelem,1)
            %Define "x" and "y"
            y = centelem(i,2);
            % parametro segundo  Starnoni
            h1=10; h2=1;
            if y>=0.5
                % solucao analitica
                u(i,1)= 11-h1*y;
                
                % calculo do gravidade
                grav(i,:)=h1*[0,1,0];
                grav_elem_escalar(i,1)=-11+h1*y;
            else
                % solucao analitica
                u(i,1)= 6.5-h2*y;
                % calculo do gravidade
                grav(i,:)=h2*[0,1,0];
                grav_elem_escalar(i,1)=-6.5+h2*y;
            end
        end
        for jj=1:size(coord,1)
            %Define "x" and "y"
            h1=10; h2=1;
            y2 = coord(jj,2);
            % parametro segundo  Starnoni
            
            if y2>=0.5
                % solucao analitica
                gravno(jj,1)= -11+h1*y2;
            else
                % solucao analitica
                gravno(jj,1)= -6.5+h2*y2;
            end
        end
        
        for j=1:size(bedge,1)+size(inedge,1)
            %Define "x" and "y"
            if j<=size(bedge,1)
                v1=bedge(j,1);
                v2=bedge(j,2);
                a=0.5*(coord(v1,:)+coord(v2,:));
                % calculo da velocidade
                IJ=coord(v1,:)-coord(v2,:);
                norma=norm(IJ);
                nij=R*IJ'/norma;
                
                ym=a(1,2);
                
                if ym>=0.5
                    h1=10;
                    V=-[-0.1*0 -0 0];
                else
                    h2=1;
                    V=-[-0.1*0 -0 0];
                end
                
                %Obtain the flow rate
                
            else
                v1=inedge(j-size(bedge,1),1);
                v2=inedge(j-size(bedge,1),2);
                a=0.5*(coord(v1,:)+coord(v2,:));
                % calculo da velocidade
                IJ=coord(v1,:)-coord(v2,:);
                norma=norm(IJ);
                nij=R*IJ'/norma;
                
                ym=a(1,2);
                if ym>=0.5
                    h1=10;
                    V=-[-0.1*0 0 0];
                else
                    h2=1;
                    V=-[-0.1*0 -0 0];
                end
            end
            F(j,1) = dot(V,nij');
            y11=a(1,2);
            % parametro segundo  Starnoni
            
            if y11>=0.5
                % solucao analitica
                gravface(j)= -11+h1*y11;
            else
                % solucao analitica
                gravface(j)= -6.5+h2*y11;
            end
        end
        
        vel=F;
        K=kmap;
        elem(:,5)=1;
    case 'starnonigrav2'
        R=[0 1 0; -1 0 0; 0 0 0];
        for i = 1:size(centelem,1)
            %Define "x" and "y"
            x = centelem(i,1);
            y = centelem(i,2);
            % solucao analitica foi calculado usando pag. 385
            % Calculo II Tom Apostol
            %  sin e cos no radianes
            u(i,1)=1+sin(x)*cos(y);
            
            % gravidade
            grav(i,:)=[-cos(x)*cos(y) sin(x)*sin(y) 0];
            grav_elem_escalar(i,1)=-1-sin(x)*cos(y);
        end
        for j=1:size(coord,1)
            %Define "x" and "y"
            x1 = coord(j,1);
            y1 = coord(j,2);
            % parametro segundo  Starnoni
            
            % solucao analitica
            gravno(j,1)=-1- sin(x1)*cos(y1);
        end
        
        for j=1:size(bedge,1)+size(inedge,1)
            %Define "x" and "y"
            if j<=size(bedge,1)
                v1=bedge(j,1);
                v2=bedge(j,2);
            else
                v1=inedge(j-size(bedge,1),1);
                v2=inedge(j-size(bedge,1),2);
            end
            a=0.5*(coord(v1,:)+coord(v2,:));
            x2=a(1,1);
            y2=a(1,2);
            gravface(j,1:3)= [-cos(x2)*cos(y2) sin(x2)*sin(y2) 0];
        end
        
        for j=1:size(bedge,1)+size(inedge,1)
            %Define "x" and "y"
            if j<=size(bedge,1)
                v1=bedge(j,1);
                v2=bedge(j,2);
                a=0.5*(coord(v1,:)+coord(v2,:));
                % calculo da velocidade
                IJ=coord(v1,:)-coord(v2,:);
            else
                v1=inedge(j-size(bedge,1),1);
                v2=inedge(j-size(bedge,1),2);
                a=0.5*(coord(v1,:)+coord(v2,:));
                % calculo da velocidade
                IJ=coord(v1,:)-coord(v2,:);
            end
            norma=norm(IJ);
            nij=R*IJ'/norma;
            % note que a velocidade analitica: velocidade_pressao +
            % velocidade_gravitacional, embora essa soma foi zero porque
            % g=-grad(p).
            V=-[cos(a(1,1))*cos(a(1,2))-0.1*sin(a(1,1))*sin(a(1,2)),...
                0.1*cos(a(1,1))*cos(a(1,2))-sin(a(1,1))*sin(a(1,2)),0]+...
                [cos(a(1,1))*cos(a(1,2))-0.1*sin(a(1,1))*sin(a(1,2)),...
                0.1*cos(a(1,1))*cos(a(1,2))-sin(a(1,1))*sin(a(1,2)),0];
            F(j,1) = dot(V,nij');
            x11=a(1,1);
            y11=a(1,2);
            % parametro segundo  Starnoni
            
            
            % solucao analitica
            gravface(j)= -1- sin(x11)*cos(y11);
            
        end
        
        vel=F;
        K=kmap;
        elem(:,5)=1;
    case 'starnonigrav3'
        R=[0 1 0; -1 0 0; 0 0 0];
        h1=10;
        h2=1;
        for i = 1:size(centelem,1)
            %Define "x" and "y"
            x = centelem(i,1);
            y = centelem(i,2);
            % parametro segundo  Starnoni
            
            if single(y)>0.5
                
                % solucao analitica
                u(i,1)= sin(x)*cos(y)+11-h1*y;
                
                % calculo do gravidade
                grav(i,:)=[-cos(x)*cos(y) h1+sin(x)*sin(y) 0];
                grav_elem_escalar(i,1)=-sin(x)*cos(y)-11+h1*y;
            else
                % solucao analitica
                u(i,1)= sin(x)*cos(y)+6.5-h2*y;
                % calculo do gravidade
                grav(i,:)=[-cos(x)*cos(y) h2+sin(x)*sin(y) 0];
                grav_elem_escalar(i,1)=-sin(x)*cos(y)-6.5+h2*y;
            end
        end
        for j=1:size(coord,1)
            %Define "x" and "y"
            x21=coord(j,1);
            y21 = coord(j,2);
            
            if single(y21)>0.5
                
                % solucao analitica
                gravno(j,1)= -sin(x21)*cos(y21)-11+h1*y21;
            else
                % solucao analitica
                gravno(j,1)= -sin(x21)*cos(y21)-6.5+h2*y21;
                
            end
            
        end
        
        for jj=1:size(bedge,1)+size(inedge,1)
            %Define "x" and "y"
            if jj<=size(bedge,1)
                v1=bedge(jj,1);
                v2=bedge(jj,2);
            else
                v1=inedge(jj-size(bedge,1),1);
                v2=inedge(jj-size(bedge,1),2);
            end
            a=0.5*(coord(v1,:)+coord(v2,:));
            x11=a(1,1);
            y11=a(1,2);
            
            if single(y11)>0.5
                
                % solucao analitica
                aaa= [-cos(x11)*cos(y11), sin(x11)*sin(y11)+h1, 0];
            else
                % solucao analitica
                aaa= [-cos(x11)*cos(y11), sin(x11)*sin(y11)+h2,0];
                
            end
            if single(y11)>0.5
                
                % solucao analitica
                gravface(jj,1)= -sin(x11)*cos(y11)-11+h1*y11;
            else
                % solucao analitica
                gravface(jj,1)= -sin(x11)*cos(y11)-6.5+h2*y11;
                
            end
            
        end
        %          for j=1:size(bedge,1)+size(inedge,1)
        %             %Define "x" and "y"
        %             if j<=size(bedge,1)
        %                 v1=bedge(j,1);
        %                 v2=bedge(j,2);
        %                 a=0.5*(coord(v1,:)+coord(v2,:));
        %                 % calculo da velocidade
        %                 IJ=coord(v1,:)-coord(v2,:);
        %             else
        %                 v1=inedge(j-size(bedge,1),1);
        %                 v2=inedge(j-size(bedge,1),2);
        %                 a=0.5*(coord(v1,:)+coord(v2,:));
        %                 % calculo da velocidade
        %                 IJ=coord(v1,:)-coord(v2,:);
        %             end
        %             y111=a(1,2);
        %             if single(y111)>0.5
        %                 V=-[cos(a(1,1))*cos(a(1,2))-0.1*(sin(a(1,1))*sin(a(1,2))+h1),...
        %                    0.1*cos(a(1,1))*cos(a(1,2))-(sin(a(1,1))*sin(a(1,2))+h1),0]+...
        %                    [cos(a(1,1))*cos(a(1,2))-0.1*(sin(a(1,1))*sin(a(1,2))+h1),...
        %                    0.1*cos(a(1,1))*cos(a(1,2))-(sin(a(1,1))*sin(a(1,2))+h1),0];
        %             else
        %                V=-[cos(a(1,1))*cos(a(1,2))-0.1*(sin(a(1,1))*sin(a(1,2))+h2),...
        %                    0.1*cos(a(1,1))*cos(a(1,2))-(sin(a(1,1))*sin(a(1,2))+h2),0]+...
        %                    [cos(a(1,1))*cos(a(1,2))-0.1*(sin(a(1,1))*sin(a(1,2))+h2),...
        %                    0.1*cos(a(1,1))*cos(a(1,2))-(sin(a(1,1))*sin(a(1,2))+h2),0];
        %             end
        %             norma=norm(IJ);
        %             nij=R*IJ'/norma;
        %             % note que a velocidade analitica: velocidade_pressao +
        %             % velocidade_gravitacional, embora essa soma sera zero porque
        %             % g=-grad(p).
        %
        %             F(j,1) = dot(V,nij');
        %         end
        %
        %         vel=F;
        K=kmap;
        elem(:,5)=1;
    case 'starnonigrav4'
        % parametro segundo  Starnoni
        h1=10;
        h2=1;
        for i = 1:size(centelem,1)
            %Define "x" and "y"
            x = centelem(i,1);
            y = centelem(i,2);
            
            if single(y)>0.5
                
                % solucao analitica
                u(i,1)= 100*sin(x)*cos(y)+11-h1*y;
                
                % calculo do gravidade
                grav(i,:)=[-100*cos(x)*cos(y) h1+100*sin(x)*sin(y) 0];
                grav_elem_escalar(i,1)=-100*sin(x)*cos(y)-11+h1*y;
            else
                % solucao analitica
                u(i,1)= 100*sin(x)*cos(y)+6.5-h2*y;
                % calculo do gravidade
                grav(i,:)=[-100*cos(x)*cos(y) h2+100*sin(x)*sin(y) 0];
                grav_elem_escalar(i,1)=-100*sin(x)*cos(y)-6.5+h2*y;
            end
        end
        for j=1:size(coord,1)
            %Define "x" and "y"
            x21=coord(j,1);
            y21 = coord(j,2);
            
            if single(y21)>0.5
                
                % solucao analitica
                gravno(j,1)= -100*sin(x21)*cos(y21)-11+h1*y21;
            else
                % solucao analitica
                gravno(j,1)= -100*sin(x21)*cos(y21)-6.5+h2*y21;
                
            end
            
        end
        
        for jj=1:size(bedge,1)+size(inedge,1)
            %Define "x" and "y"
            if jj<=size(bedge,1)
                v1=bedge(jj,1);
                v2=bedge(jj,2);
            else
                v1=inedge(jj-size(bedge,1),1);
                v2=inedge(jj-size(bedge,1),2);
            end
            a=0.5*(coord(v1,:)+coord(v2,:));
            x11=a(1,1);
            y11=a(1,2);
            
            if single(y11)>0.5
                bbb= [-100*cos(x11)*cos(y11), 100*sin(x11)*sin(y11)+h1];
            else
                bbb= [-100*cos(x11)*cos(y11), 100*sin(x11)*sin(y11)+h2];
                
            end
            if single(y11)>0.5
                
                % solucao analitica
                gravface(jj,1)= -100*sin(x11)*cos(y11)-11+h1*y11;
            else
                % solucao analitica
                gravface(jj,1)= -100*sin(x11)*cos(y11)-6.5+h2*y11;
                
            end
        end
        
        K=kmap;
        elem(:,5)=1;
    case 'zhangkobaise'
        R=[0 1 0; -1 0 0; 0 0 0];
        for i = 1:size(centelem,1)
            
            %Define "x" and "y"
            x = centelem(i,1);
            y = centelem(i,2);
            teta=calculoteta(x,y);
            alfa=0.5645;
            r=sqrt(x^2+y^2);
            coef=alfa*(alfa-1)*(x^2+y^2)^(0.5*alfa-2);
            
            if elem(i,5)==1
                a1=1.0;
                b1=-12.0414;
                
                u(i,1)= 10+(r^alfa)*(a1*cos(alfa*teta)+b1*sin(alfa*teta));
                kmap(1,:) = [1 1 0 0 1];
                fonte(i,1) =0;
            elseif elem(i,5)==2
                a2=-4.8591;
                b2=-6.0699;
                
                u(i,1)= 10+(r^alfa)*(a2*cos(alfa*teta)+b2*sin(alfa*teta));
                kmap(2,:) = [2 10 0 0 10];
                %==========================================================
                sum1= (x^2)*(a2*cos(alfa*teta)+b2*sin(alfa*teta));
                sum2=(x*y)*(-2*b2*cos(alfa*teta)+2*a2*sin(alfa*teta));
                sum3=(y^2)*(-a2*cos(alfa*teta)-b2*sin(alfa*teta));
                sumtotal=sum1+sum2+sum3;
                fonte(i,1) =(10-10)*coef*sumtotal*elemarea(i,1);
            else
                a3=-0.9664;
                b3=-0.2837;
                
                u(i,1)= 10+(r^alfa)*(a3*cos(alfa*teta)+b3*sin(alfa*teta));
                kmap(3,:) = [3 100 0 0 100];
                fonte(i,1) =0;
            end
            
            
        end  %End of FOR
        K=kmap;
        
        %% velocities
        for ianalit = 1:size(bedge,1)
            lef=bedge(ianalit,3);
            IJ=coord(bedge(ianalit,2),:)-coord(bedge(ianalit,1),:);
            norma=norm(IJ);
            nij=R*IJ'/norma;
            
            auxpoint=(coord(bedge(ianalit,2),:)+coord(bedge(ianalit,1),:))*0.5;
            x=auxpoint(1,1);
            y=auxpoint(1,2);
            
            teta1=calculoteta(x,y);
            alfa=0.5645;
            r=sqrt(x^2+y^2);
            
            if elem(lef,5)==1
                a1=1.0;
                b1=-12.0414;
                k11=1;
                k22=1;
                dtetadx=calculotetadx(x,y);
                part1=((alfa*(r^alfa))/r^2)*x*(a1*cos(alfa*teta1)+b1*sin(alfa*teta1));
                part2=alfa*(r^alfa)*(-a1*sin(alfa*teta1)+b1*cos(alfa*teta1))*dtetadx;
                dtetady=calculotetady(x,y);
                part3=((alfa*(r^alfa))/r^2)*y*(a1*cos(alfa*teta1)+b1*sin(alfa*teta1));
                part4=alfa*(r^alfa)*(-a1*sin(alfa*teta1)+b1*cos(alfa*teta1))*dtetady;
                
                V=-[k11*(part1+part2) k22*(part3+part4) 0];
            elseif elem(lef,5)==2
                a2=-4.8591;
                b2=-6.0699;
                k11=10;
                k22=10;
                % dtheta/dx
                dtetadx=calculotetadx(x,y);
                %==========================================================
                var1=((alfa*(r^alfa))/r^2)*(a2*cos(alfa*teta1)+b2*sin(alfa*teta1));
                var2=alfa*(r^alfa)*(-a2*sin(alfa*teta1)+b2*cos(alfa*teta1));
                %
                part1=var1*x;
                part2=var2*dtetadx;
                % dtheta/dy
                dtetady=calculotetady(x,y);
                part3=var1*y;
                part4=var2*dtetady;
                V=-[k11*(part1+part2) k22*(part3+part4) 0];
            else
                a3=-0.9664;
                b3=-0.2837;
                k11=100;
                k22=100;
                dtetadx=calculotetadx(x,y);
                part1=((alfa*(r^alfa))/r^2)*x*(a3*cos(alfa*teta1)+b3*sin(alfa*teta1));
                part2=alfa*(r^alfa)*(-a3*sin(alfa*teta1)+b3*cos(alfa*teta1))*dtetadx;
                dtetady=calculotetady(x,y);
                part3=((alfa*(r^alfa))/r^2)*y*(a3*cos(alfa*teta1)+b3*sin(alfa*teta1));
                part4=alfa*(r^alfa)*(-a3*sin(alfa*teta1)+b3*cos(alfa*teta1))*dtetady;
                V=-[k11*(part1+part2) k22*(part3+part4) 0];
            end
            %Obtain the flow rate
            F(ianalit,1) = dot(V,nij');
        end  %End of FOR (boundary edges)
        
        for ianalit=1:size(inedge,1)
            lef=inedge(ianalit,3);
            
            IJ=coord(inedge(ianalit,2),:)-coord(inedge(ianalit,1),:);
            norma=norm(IJ);
            nij=R*IJ'/norma;
            auxpoint=(coord(inedge(ianalit,2),:)+coord(inedge(ianalit,1),:))*0.5;
            x=auxpoint(1,1);
            y=auxpoint(1,2);
            
            
            teta1=calculoteta(x,y);
            alfa=0.5645;
            r=sqrt(x^2+y^2);
            
            if elem(lef,5)==1
                a1=1.0;
                b1=-12.0414;
                k11=1;
                k22=1;
                dtetadx=calculotetadx(x,y);
                part1=((alfa*(r^alfa))/r^2)*x*(a1*cos(alfa*teta1)+b1*sin(alfa*teta1));
                part2=alfa*(r^alfa)*(-a1*sin(alfa*teta1)+b1*cos(alfa*teta1))*dtetadx;
                dtetady=calculotetady(x,y);
                part3=((alfa*(r^alfa))/r^2)*y*(a1*cos(alfa*teta1)+b1*sin(alfa*teta1));
                part4=alfa*(r^alfa)*(-a1*sin(alfa*teta1)+b1*cos(alfa*teta1))*dtetady;
                
                V=-[k11*(part1+part2) k22*(part3+part4) 0];
            elseif elem(lef,5)==2
                a2=-4.8591;
                b2=-6.0699;
                k11=10;
                k22=10;
                dtetadx=calculotetadx(x,y);
                part1=((alfa*(r^alfa))/r^2)*x*(a2*cos(alfa*teta1)+b2*sin(alfa*teta1));
                part2=alfa*(r^alfa)*(-a2*sin(alfa*teta1)+b2*cos(alfa*teta1))*dtetadx;
                dtetady=calculotetady(x,y);
                part3=((alfa*(r^alfa))/r^2)*y*(a2*cos(alfa*teta1)+b2*sin(alfa*teta1));
                part4=alfa*(r^alfa)*(-a2*sin(alfa*teta1)+b2*cos(alfa*teta1))*dtetady;
                V=-[k11*(part1+part2) k22*(part3+part4) 0];
            else
                a3=-0.9664;
                b3=-0.2837;
                k11=100;
                k22=100;
                dtetadx=calculotetadx(x,y);
                part1=((alfa*(r^alfa))/r^2)*x*(a3*cos(alfa*teta1)+b3*sin(alfa*teta1));
                part2=alfa*(r^alfa)*(-a3*sin(alfa*teta1)+b3*cos(alfa*teta1))*dtetadx;
                dtetady=calculotetady(x,y);
                part3=((alfa*(r^alfa))/r^2)*y*(a3*cos(alfa*teta1)+b3*sin(alfa*teta1));
                part4=alfa*(r^alfa)*(-a3*sin(alfa*teta1)+b3*cos(alfa*teta1))*dtetady;
                V=-[k11*(part1+part2) k22*(part3+part4) 0];
            end
            
            F(ianalit+size(bedge,1),1) = dot(V,nij');
        end
        
        vel=F;
    case 'zhangkobaise2'
        R=[0 1 0; -1 0 0; 0 0 0];
        for i = 1:size(centelem,1)
            
            %Define "x" and "y"
            x = centelem(i,1);
            y = centelem(i,2);
            teta=calculoteta(x,y);
            alfa=0.6142;
            r=sqrt(x^2+y^2);
            coef=alfa*(alfa-1)*(x^2+y^2)^(0.5*alfa-2);
            if elem(i,5)==1
                a1=1.0;
                b1=-1.0546;
                
                u(i,1)= 10+(r^alfa)*(a1*cos(alfa*teta)+b1*sin(alfa*teta));
                kmap(1,:) = [1 1 0 0 1];
                fonte(i,1) =0;
            elseif elem(i,5)==2
                a2=-0.4275;
                b2=0.2142;
                
                u(i,1)= 10+(r^alfa)*(a2*cos(alfa*teta)+b2*sin(alfa*teta));
                kmap(2,:) = [2 10 0 0 1000];
                %==========================================================
                sum1= (x^2)*(a2*cos(alfa*teta)+b2*sin(alfa*teta));
                sum2=(x*y)*(-2*b2*cos(alfa*teta)+2*a2*sin(alfa*teta));
                sum3=(y^2)*(-a2*cos(alfa*teta)-b2*sin(alfa*teta));
                sumtotal=sum1+sum2+sum3;
                fonte(i,1) =(10-1000)*coef*sumtotal*elemarea(i,1);
            else
                a3=-0.7604;
                b3=-0.6495;
                
                u(i,1)= 10+(r^alfa)*(a3*cos(alfa*teta)+b3*sin(alfa*teta));
                kmap(3,:) = [3 100 0 0 100];
                fonte(i,1) =0;
            end
            
            
        end  %End of FOR
        K=kmap;
        
        %% velocities
        for ianalit = 1:size(bedge,1)
            lef=bedge(ianalit,3);
            IJ=coord(bedge(ianalit,2),:)-coord(bedge(ianalit,1),:);
            norma=norm(IJ);
            nij=R*IJ'/norma;
            
            auxpoint=(coord(bedge(ianalit,2),:)+coord(bedge(ianalit,1),:))*0.5;
            x=auxpoint(1,1);
            y=auxpoint(1,2);
            
            teta1=calculoteta(x,y);
            alfa=0.6142;
            r=sqrt(x^2+y^2);
            
            if elem(lef,5)==1
                a1=1.0;
                b1=-1.0546;
                k11=1;
                k22=1;
                dtetadx=calculotetadx(x,y);
                part1=((alfa*(r^alfa))/r^2)*x*(a1*cos(alfa*teta1)+b1*sin(alfa*teta1));
                part2=alfa*(r^alfa)*(-a1*sin(alfa*teta1)+b1*cos(alfa*teta1))*dtetadx;
                dtetady=calculotetady(x,y);
                part3=((alfa*(r^alfa))/r^2)*y*(a1*cos(alfa*teta1)+b1*sin(alfa*teta1));
                part4=alfa*(r^alfa)*(-a1*sin(alfa*teta1)+b1*cos(alfa*teta1))*dtetady;
                
                V=-[k11*(part1+part2) k22*(part3+part4) 0];
            elseif elem(lef,5)==2
                a2=-0.4275;
                b2=0.2142;
                k11=10;
                k22=1000;
                % dtheta/dx
                dtetadx=calculotetadx(x,y);
                %==========================================================
                var1=((alfa*(r^alfa))/r^2)*(a2*cos(alfa*teta1)+b2*sin(alfa*teta1));
                var2=alfa*(r^alfa)*(-a2*sin(alfa*teta1)+b2*cos(alfa*teta1));
                %
                part1=var1*x;
                part2=var2*dtetadx;
                % dtheta/dy
                dtetady=calculotetady(x,y);
                part3=var1*y;
                part4=var2*dtetady;
                V=-[k11*(part1+part2) k22*(part3+part4) 0];
            else
                a3=-0.7604;
                b3=-0.6495;
                k11=100;
                k22=100;
                dtetadx=calculotetadx(x,y);
                part1=((alfa*(r^alfa))/r^2)*x*(a3*cos(alfa*teta1)+b3*sin(alfa*teta1));
                part2=alfa*(r^alfa)*(-a3*sin(alfa*teta1)+b3*cos(alfa*teta1))*dtetadx;
                dtetady=calculotetady(x,y);
                part3=((alfa*(r^alfa))/r^2)*y*(a3*cos(alfa*teta1)+b3*sin(alfa*teta1));
                part4=alfa*(r^alfa)*(-a3*sin(alfa*teta1)+b3*cos(alfa*teta1))*dtetady;
                V=-[k11*(part1+part2) k22*(part3+part4) 0];
            end
            %Obtain the flow rate
            F(ianalit,1) = dot(V,nij');
        end  %End of FOR (boundary edges)
        
        for ianalit=1:size(inedge,1)
            lef=inedge(ianalit,3);
            
            IJ=coord(inedge(ianalit,2),:)-coord(inedge(ianalit,1),:);
            norma=norm(IJ);
            nij=R*IJ'/norma;
            auxpoint=(coord(inedge(ianalit,2),:)+coord(inedge(ianalit,1),:))*0.5;
            x=auxpoint(1,1);
            y=auxpoint(1,2);
            
            
            teta1=calculoteta(x,y);
            alfa=0.6142;
            r=sqrt(x^2+y^2);
            
            if elem(lef,5)==1
                a1=1.0;
                b1=-1.0546;
                k11=1;
                k22=1;
                dtetadx=calculotetadx(x,y);
                part1=((alfa*(r^alfa))/r^2)*x*(a1*cos(alfa*teta1)+b1*sin(alfa*teta1));
                part2=alfa*(r^alfa)*(-a1*sin(alfa*teta1)+b1*cos(alfa*teta1))*dtetadx;
                dtetady=calculotetady(x,y);
                part3=((alfa*(r^alfa))/r^2)*y*(a1*cos(alfa*teta1)+b1*sin(alfa*teta1));
                part4=alfa*(r^alfa)*(-a1*sin(alfa*teta1)+b1*cos(alfa*teta1))*dtetady;
                
                V=-[k11*(part1+part2) k22*(part3+part4) 0];
            elseif elem(lef,5)==2
                a2=-0.4275;
                b2=0.2142;
                k11=10;
                k22=1000;
                dtetadx=calculotetadx(x,y);
                part1=((alfa*(r^alfa))/r^2)*x*(a2*cos(alfa*teta1)+b2*sin(alfa*teta1));
                part2=alfa*(r^alfa)*(-a2*sin(alfa*teta1)+b2*cos(alfa*teta1))*dtetadx;
                dtetady=calculotetady(x,y);
                part3=((alfa*(r^alfa))/r^2)*y*(a2*cos(alfa*teta1)+b2*sin(alfa*teta1));
                part4=alfa*(r^alfa)*(-a2*sin(alfa*teta1)+b2*cos(alfa*teta1))*dtetady;
                V=-[k11*(part1+part2) k22*(part3+part4) 0];
            else
                a3=-0.7604;
                b3=-0.6495;
                k11=100;
                k22=100;
                dtetadx=calculotetadx(x,y);
                part1=((alfa*(r^alfa))/r^2)*x*(a3*cos(alfa*teta1)+b3*sin(alfa*teta1));
                part2=alfa*(r^alfa)*(-a3*sin(alfa*teta1)+b3*cos(alfa*teta1))*dtetadx;
                dtetady=calculotetady(x,y);
                part3=((alfa*(r^alfa))/r^2)*y*(a3*cos(alfa*teta1)+b3*sin(alfa*teta1));
                part4=alfa*(r^alfa)*(-a3*sin(alfa*teta1)+b3*cos(alfa*teta1))*dtetady;
                V=-[k11*(part1+part2) k22*(part3+part4) 0];
            end
            
            F(ianalit+size(bedge,1),1) = dot(V,nij');
        end
        
        vel=F;
        
        
    case 'zhangkobaise3'
        R=[0 1 0; -1 0 0; 0 0 0];
        for i = 1:size(centelem,1)
            
            %Define "x" and "y"
            x = centelem(i,1);
            y = centelem(i,2);
            teta=calculoteta(x,y);
            alfa=0.8866;
            r=sqrt(x^2+y^2);
            coef=alfa*(alfa-1)*(x^2+y^2)^(0.5*alfa-2);
            if elem(i,5)==1
                a1=1.0;
                b1=-0.3706;
                
                u(i,1)= 10+(r^alfa)*(a1*cos(alfa*teta)+b1*sin(alfa*teta));
                kmap(1,:) = [1 1 0 0 1];
                fonte(i,1) =0;
            elseif elem(i,5)==2
                a2=-0.0144;
                b2=0.0022;
                
                u(i,1)= 10+(r^alfa)*(a2*cos(alfa*teta)+b2*sin(alfa*teta));
                kmap(2,:) = [2 10 0 0 100000];
                %=========================================================
                sum1= (x^2)*(a2*cos(alfa*teta)+b2*sin(alfa*teta));
                sum2=(x*y)*(-2*b2*cos(alfa*teta)+2*a2*sin(alfa*teta));
                sum3=(y^2)*(-a2*cos(alfa*teta)-b2*sin(alfa*teta));
                sumtotal=sum1+sum2+sum3;
                fonte(i,1) =(10-100000)*coef*sumtotal*elemarea(i,1);
            else
                a3=0.7544;
                b3=-0.6564;
                
                u(i,1)= 10+(r^alfa)*(a3*cos(alfa*teta)+b3*sin(alfa*teta));
                kmap(3,:) = [3 100 0 0 100];
                fonte(i,1) =0;
            end
        end  %End of FOR
        K=kmap;
        %% velocities
        for ianalit = 1:size(bedge,1)
            lef=bedge(ianalit,3);
            IJ=coord(bedge(ianalit,2),:)-coord(bedge(ianalit,1),:);
            norma=norm(IJ);
            nij=R*IJ'/norma;
            
            auxpoint=(coord(bedge(ianalit,2),:)+coord(bedge(ianalit,1),:))*0.5;
            x=auxpoint(1,1);
            y=auxpoint(1,2);
            
            teta1=calculoteta(x,y);
            alfa=0.8866;
            r=sqrt(x^2+y^2);
            
            if elem(lef,5)==1
                a1=1.0;
                b1=-0.3706;
                k11=1;
                k22=1;
                dtetadx=calculotetadx(x,y);
                part1=((alfa*(r^alfa))/r^2)*x*(a1*cos(alfa*teta1)+b1*sin(alfa*teta1));
                part2=alfa*(r^alfa)*(-a1*sin(alfa*teta1)+b1*cos(alfa*teta1))*dtetadx;
                dtetady=calculotetady(x,y);
                part3=((alfa*(r^alfa))/r^2)*y*(a1*cos(alfa*teta1)+b1*sin(alfa*teta1));
                part4=alfa*(r^alfa)*(-a1*sin(alfa*teta1)+b1*cos(alfa*teta1))*dtetady;
                
                V=-[k11*(part1+part2) k22*(part3+part4) 0];
            elseif elem(lef,5)==2
                a2=-0.0144;
                b2=0.0022;
                k11=10;
                k22=100000;
                % dtheta/dx
                dtetadx=calculotetadx(x,y);
                %==========================================================
                var1=((alfa*(r^alfa))/r^2)*(a2*cos(alfa*teta1)+b2*sin(alfa*teta1));
                var2=alfa*(r^alfa)*(-a2*sin(alfa*teta1)+b2*cos(alfa*teta1));
                %
                part1=var1*x;
                part2=var2*dtetadx;
                % dtheta/dy
                dtetady=calculotetady(x,y);
                part3=var1*y;
                part4=var2*dtetady;
                V=-[k11*(part1+part2) k22*(part3+part4) 0];
            else
                a3=0.7544;
                b3=-0.6564;
                k11=100;
                k22=100;
                dtetadx=calculotetadx(x,y);
                part1=((alfa*(r^alfa))/r^2)*x*(a3*cos(alfa*teta1)+b3*sin(alfa*teta1));
                part2=alfa*(r^alfa)*(-a3*sin(alfa*teta1)+b3*cos(alfa*teta1))*dtetadx;
                dtetady=calculotetady(x,y);
                part3=((alfa*(r^alfa))/r^2)*y*(a3*cos(alfa*teta1)+b3*sin(alfa*teta1));
                part4=alfa*(r^alfa)*(-a3*sin(alfa*teta1)+b3*cos(alfa*teta1))*dtetady;
                V=-[k11*(part1+part2) k22*(part3+part4) 0];
            end
            %Obtain the flow rate
            F(ianalit,1) = dot(V,nij');
        end  %End of FOR (boundary edges)
        
        for ianalit=1:size(inedge,1)
            lef=inedge(ianalit,3);
            
            IJ=coord(inedge(ianalit,2),:)-coord(inedge(ianalit,1),:);
            norma=norm(IJ);
            nij=R*IJ'/norma;
            auxpoint=(coord(inedge(ianalit,2),:)+coord(inedge(ianalit,1),:))*0.5;
            x=auxpoint(1,1);
            y=auxpoint(1,2);
            
            
            teta1=calculoteta(x,y);
            alfa=0.8866;
            r=sqrt(x^2+y^2);
            
            if elem(lef,5)==1
                a1=1.0;
                b1=-0.3706;
                k11=1;
                k22=1;
                dtetadx=calculotetadx(x,y);
                part1=((alfa*(r^alfa))/r^2)*x*(a1*cos(alfa*teta1)+b1*sin(alfa*teta1));
                part2=alfa*(r^alfa)*(-a1*sin(alfa*teta1)+b1*cos(alfa*teta1))*dtetadx;
                dtetady=calculotetady(x,y);
                part3=((alfa*(r^alfa))/r^2)*y*(a1*cos(alfa*teta1)+b1*sin(alfa*teta1));
                part4=alfa*(r^alfa)*(-a1*sin(alfa*teta1)+b1*cos(alfa*teta1))*dtetady;
                
                V=-[k11*(part1+part2) k22*(part3+part4) 0];
            elseif elem(lef,5)==2
                a2=-0.0144;
                b2=0.0022;
                k11=10;
                k22=100000;
                dtetadx=calculotetadx(x,y);
                part1=((alfa*(r^alfa))/r^2)*x*(a2*cos(alfa*teta1)+b2*sin(alfa*teta1));
                part2=alfa*(r^alfa)*(-a2*sin(alfa*teta1)+b2*cos(alfa*teta1))*dtetadx;
                dtetady=calculotetady(x,y);
                part3=((alfa*(r^alfa))/r^2)*y*(a2*cos(alfa*teta1)+b2*sin(alfa*teta1));
                part4=alfa*(r^alfa)*(-a2*sin(alfa*teta1)+b2*cos(alfa*teta1))*dtetady;
                V=-[k11*(part1+part2) k22*(part3+part4) 0];
            else
                a3=0.7544;
                b3=-0.6564;
                k11=100;
                k22=100;
                dtetadx=calculotetadx(x,y);
                part1=((alfa*(r^alfa))/r^2)*x*(a3*cos(alfa*teta1)+b3*sin(alfa*teta1));
                part2=alfa*(r^alfa)*(-a3*sin(alfa*teta1)+b3*cos(alfa*teta1))*dtetadx;
                dtetady=calculotetady(x,y);
                part3=((alfa*(r^alfa))/r^2)*y*(a3*cos(alfa*teta1)+b3*sin(alfa*teta1));
                part4=alfa*(r^alfa)*(-a3*sin(alfa*teta1)+b3*cos(alfa*teta1))*dtetady;
                V=-[k11*(part1+part2) k22*(part3+part4) 0];
            end
            
            F(ianalit+size(bedge,1),1) = dot(V,nij');
        end
        
        vel=F;
        
    case 'zigzagfract'
        for i = 1:size(centelem,1)
            
            %Define "x" and "y"
            x = centelem(i,1);
            y = centelem(i,2);
            if x<0.292
                if (0.15<x && x<0.292) && (0.48<y && y<0.52)
                    elem(i,5)=4;
                    
                    kmap(4,:) = [4 0.75025e-3 0.43258e-3 0.43258e-3 0.25025e-3];
                    %fonte(i,1)=100;
                else
                    elem(i,5)=1;
                    
                    kmap(1,:) = [1 0.75025e-3 0.43258e-3 0.43258e-3 0.25025e-3];
                    %fonte(i,1)=0;
                end
            elseif 0.292<x && x<0.711
                if (0.292<x && x<0.711) && (0.48<y && y<0.52)
                    elem(i,5)=5;
                    
                    kmap(5,:) = [5 0.25025e-3 -0.43258e-3 -0.43258e-3 0.75025e-3];
                    %fonte(i,1)=100;
                else
                    elem(i,5)=2;
                    
                    kmap(2,:) = [2 0.25025e-3 -0.43258e-3 -0.43258e-3 0.75025e-3];
                    %fonte(i,1)=0;
                end
            else
                
                if ((0.711<x && x<0.85) && (0.48<y && y<0.52))
                    elem(i,5)=6;
                    
                    kmap(6,:) = [6 0.25025e-3 0.43258e-3 0.43258e-3 0.75025e-3];
                    %fonte(i,1)=100;
                else
                    elem(i,5)=3;
                    
                    kmap(3,:) = [3 0.25025e-3 0.43258e-3 0.43258e-3 0.75025e-3];
                    %fonte(i,1)=0;
                end
            end
            
        end  %End of FOR
        K=kmap;
        
    case 'herbinhubert'
        
        %Initialize "kmap"
        kmap = zeros(size(centelem,1),5);
        for i = 1:size(centelem,1)
            elem(i,5)=1;
            %Define "x" and "y"
            x = centelem(i,1);
            y = centelem(i,2);
            %Definition of permeability components
            k(1,1) = 1.5;
            k(1,2) = 0.5;
            k(2,1) = 0.5;
            k(2,2) = 1.5;
            
            %Build "kmap"
            kmap(i,:) = [i k(1,1) k(1,2) k(2,1) k(2,2)];
            u(i,1)= 16*x*(1-x)*y*(1-y);
            fonte(i,1) =1.5*32*(y*(1-y)+x*(1-x))-16*(1-2*x)*(1-2*y);
        end  %End of FOR
        K=kmap;
        
    case 'herbin'
        R=[0 -1 0; 1 0 0; 0 0 0];
        %Initialize a parameters
        alfa = 1000;
        %Initialize "kmap"
        kmap = zeros(size(centelem,1),5);
        for i = 1:size(centelem,1)
            elem(i,5)=i;
            %Define "x" and "y"
            x = centelem(i,1);
            y = centelem(i,2);
            %Definition of permeability components
            kaux(1,1) = (alfa*(x^2) + (y^2));
            kaux(1,2) = (alfa - 1)*(x*y);
            kaux(2,1) = (alfa - 1)*(x*y);
            kaux(2,2) = (alfa*(y^2) + (x^2));
            
            k = (1/((x^2) + (y^2))).*kaux;
            
            %Build "kmap"
            kmap(i,:) = [i k(1,1) k(1,2) k(2,1) k(2,2)];
            u(i,1)= sin(pi*x)*sin(pi*y);
            fonte(i,1) =(pi/((centelem(i,1)^2) + ...
                (centelem(i,2)^2)))*(-(alfa - 1)*centelem(i,1)*...
                cos(pi*centelem(i,1))*(2*pi*centelem(i,2)*...
                cos(pi*centelem(i,2)) + sin(pi*centelem(i,2))) + ...
                sin(pi*centelem(i,1))*(-(alfa - 1)*centelem(i,2)*...
                cos(pi*centelem(i,2)) + (1 + alfa)*pi*...
                ((centelem(i,1)^2) + (centelem(i,2)^2))*...
                sin(pi*centelem(i,2))));
            
        end  %End of FOR
        K=kmap;
        %Initialize "velanalit"
        F = zeros(size(bedge,1) + size(inedge,1),1);
        %Calculate the analitical velocity
        %Boundary edges
        for ianalit = 1:size(bedge,1)
            
            %Obtain the vector velocity (-KnablaP)
            %     V = [(-pi/((overedgecoord(ianal,1)^2) + (overedgecoord(ianal,2)^2)))*...
            %         ((delta - 1)*overedgecoord(ianal,1)*overedgecoord(ianal,2)*...
            %         cos(pi*overedgecoord(ianal,2))*sin(pi*overedgecoord(ianal,1)) + ...
            %         (delta*(overedgecoord(ianal,1)^2) + (overedgecoord(ianal,2)^2))*...
            %         cos(pi*overedgecoord(ianal,1))*sin(pi*overedgecoord(ianal,2)));
            %         (-pi/((overedgecoord(ianal,1)^2) + (overedgecoord(ianal,2)^2)))*...
            %         ((delta - 1)*overedgecoord(ianal,1)*overedgecoord(ianal,2)*...
            %         cos(pi*overedgecoord(ianal,1))*sin(pi*overedgecoord(ianal,2)) + ...
            %         (delta*(overedgecoord(ianal,2)^2) + (overedgecoord(ianal,1)^2))*...
            %         cos(pi*overedgecoord(ianal,2))*sin(pi*overedgecoord(ianal,1)));
            %         0];
            IJ=coord(bedge(ianalit,2),:)-coord(bedge(ianalit,1),:);
            norma=norm(IJ);
            nij=R*IJ'/norma;
            
            auxpoint=(coord(bedge(ianalit,2),:)+coord(bedge(ianalit,1),:))*0.5;
            x=auxpoint(1,1);
            y=auxpoint(1,2);
            
            dpdx=pi*cos(pi*x)*sin(pi*y);
            dpdy=pi*sin(pi*x)*cos(pi*y);
            k11=(alfa*x^2 + y^2)/(x^2+y^2);
            k12= ((alfa-1)*x*y)/(x^2+y^2);
            k22=(x^2+alfa*y^2)/(x^2+y^2);
            
            V=-[k11*dpdx+k12*dpdy k12*dpdx+k22*dpdy 0];
            %Obtain the flow rate
            F(ianalit) = dot(V,nij');
        end  %End of FOR (boundary edges)
        
        %Internal edges
        for ianalit = 1:size(inedge,1)
            IJ=coord(inedge(ianalit,2),:)-coord(inedge(ianalit,1),:);
            norma=norm(IJ);
            nij=R*IJ'/norma;
            auxpoint=(coord(inedge(ianalit,2),:)+coord(inedge(ianalit,1),:))*0.5;
            x=auxpoint(1,1);
            y=auxpoint(1,2);
            
            dpdx=pi*cos(pi*x)*sin(pi*y);
            dpdy=pi*sin(pi*x)*cos(pi*y);
            k11=(alfa*x^2 + y^2)/(x^2+y^2);
            k12= ((alfa-1)*x*y)/(x^2+y^2);
            k22=(x^2+alfa*y^2)/(x^2+y^2);
            
            V=-[k11*dpdx+k12*dpdy k12*dpdx+k22*dpdy 0];
            %Obtain the flow rate
            F(size(bedge,1) + ianalit) = ...
                dot(V,nij');
        end
        vel=F;
    case 'gaowu6'
        fonte=0;
        elem(:,5)=1;
        u=0;
        %Initialize "R":
        R = zeros(2);
        % problema Gao e Wu 2013
        k = [1 0; 0 1e-3];
        % problema TEREKHOV
        %k = [50 0; 0 1];
        %Fill "R"
        R(1,1) = cosd(67.5);
        R(1,2) = sind(67.5);
        R(2,1) = -R(1,2);
        R(2,2) = R(1,1);
        %Fill "k" turning the tensor
        A=inv(R);
        k = A*k*R;
        %Buld "kmap" again
        kmap = [1 k(1,1) k(1,2) k(2,1) k(2,2)];
        K=kmap;
    case 'benchmar5_6'
        %% problem 5.6 CASO 2
        u=0;
        for ielem=1:size(elem,1)
            if (0<centelem(ielem,1) && centelem(ielem,1)<0.5) && ...
                    (0<centelem(ielem,2) && centelem(ielem,2)<0.5)
                elem(ielem,5)=ielem;
                theta=pi/6;
                k=[1000 0;0 1];
                karot=[ cos(theta) -sin(theta); sin(theta) cos(theta)]*k*...
                    [cos(theta) sin(theta); -sin(theta) cos(theta)];
                kmap(ielem,1:5)=[ielem karot(1,1) karot(1,2) karot(2,1) karot(2,2)];
            elseif (0.5<centelem(ielem,1) && centelem(ielem,1)<1) && ...
                    (0<centelem(ielem,2) && centelem(ielem,2)<0.5)
                elem(ielem,5)=ielem;
                theta=-pi/6;
                k=[10 0;0 1];
                karot=[ cos(theta) -sin(theta); sin(theta) cos(theta)]*k*...
                    [cos(theta) sin(theta); -sin(theta) cos(theta)];
                kmap(ielem,1:5)=[ielem karot(1,1) karot(1,2) karot(2,1) karot(2,2)];
            elseif (0.5<centelem(ielem,1) && centelem(ielem,1)<1) && ...
                    (0.5<centelem(ielem,2) && centelem(ielem,2)<1)
                elem(ielem,5)=ielem;
                theta=pi/6;
                k=[1000 0;0 1];
                karot=[ cos(theta) -sin(theta); sin(theta) cos(theta)]*k*...
                    [cos(theta) sin(theta); -sin(theta) cos(theta)];
                kmap(ielem,1:5)=[ielem karot(1,1) karot(1,2) karot(2,1) karot(2,2)];
            elseif (0<centelem(ielem,1) && centelem(ielem,1)<0.5) && ...
                    (0.5<centelem(ielem,2) && centelem(ielem,2)<1)
                elem(ielem,5)=ielem;
                theta=-pi/6;
                k=[10 0;0 1];
                karot=[ cos(theta) -sin(theta); sin(theta) cos(theta)]*k*...
                    [cos(theta) sin(theta); -sin(theta) cos(theta)];
                kmap(ielem,1:5)=[ielem karot(1,1) karot(1,2) karot(2,1) karot(2,2)];
            end
            if (7/18<centelem(ielem,1) || 7/18==centelem(ielem,1))&& ...
                    (centelem(ielem,1)<11/18 || 11/18==centelem(ielem,1))...
                    && (7/18<centelem(ielem,2) || centelem(ielem,2)==7/18) && ...
                    (centelem(ielem,2)<11/18 || centelem(ielem,2)==11/18)
                fonte(ielem,1)=81/4;
            else
                fonte(ielem,1)=0;
            end
            K=kmap;
            
        end
        fonte=fonte.*elemarea;
    case 'benchmar5_7' %
        %% problem 5.6 CASO 1
        u=0;
        k=[1000 0;0 1];
        for ielem=1:size(elem,1)
            if (0<centelem(ielem,1) && centelem(ielem,1)<0.5) && ...
                    (0<centelem(ielem,2) && centelem(ielem,2)<0.5)
                elem(ielem,5)=ielem;
                theta=pi/6;
                
                karot=[ cos(theta) -sin(theta); sin(theta) cos(theta)]*k*...
                    [cos(theta) sin(theta); -sin(theta) cos(theta)];
                kmap(ielem,1:5)=[ielem karot(1,1) karot(1,2) karot(2,1) karot(2,2)];
            elseif (0.5<centelem(ielem,1) && centelem(ielem,1)<1) && ...
                    (0<centelem(ielem,2) && centelem(ielem,2)<0.5)
                elem(ielem,5)=ielem;
                theta=-pi/6;
                
                karot=[ cos(theta) -sin(theta); sin(theta) cos(theta)]*k*...
                    [cos(theta) sin(theta); -sin(theta) cos(theta)];
                kmap(ielem,1:5)=[ielem karot(1,1) karot(1,2) karot(2,1) karot(2,2)];
            elseif (0.5<centelem(ielem,1) && centelem(ielem,1)<1) && ...
                    (0.5<centelem(ielem,2) && centelem(ielem,2)<1)
                elem(ielem,5)=ielem;
                theta=pi/6;
                
                karot=[ cos(theta) -sin(theta); sin(theta) cos(theta)]*k*...
                    [cos(theta) sin(theta); -sin(theta) cos(theta)];
                kmap(ielem,1:5)=[ielem karot(1,1) karot(1,2) karot(2,1) karot(2,2)];
            elseif (0<centelem(ielem,1) && centelem(ielem,1)<0.5) && ...
                    (0.5<centelem(ielem,2) && centelem(ielem,2)<1)
                elem(ielem,5)=ielem;
                theta=-pi/6;
                
                karot=[ cos(theta) -sin(theta); sin(theta) cos(theta)]*k*...
                    [cos(theta) sin(theta); -sin(theta) cos(theta)];
                kmap(ielem,1:5)=[ielem karot(1,1) karot(1,2) karot(2,1) karot(2,2)];
            end
            if (7/18<centelem(ielem,1) || 7/18==centelem(ielem,1))&& ...
                    (centelem(ielem,1)<11/18 || 11/18==centelem(ielem,1))...
                    && (7/18<centelem(ielem,2) || centelem(ielem,2)==7/18) && ...
                    (centelem(ielem,2)<11/18 || centelem(ielem,2)==11/18)
                fonte(ielem,1)=81/4;
            else
                fonte(ielem,1)=0;
            end
        end
        K=kmap;
        fonte=fonte.*elemarea;
    case 'edwards'
        %% This problem is adapted of (Rogers and Edwards 1998)
        a1=1;
        f=((4*a1)/(((1/50)-2)*(1/10)+1));
        b2=((1/10)-1)*f;
        c2=f;
        d2=-c2*(1/10);
        c1=(1/50)*(1/10)*f;
        d1=d2;
        % calculo da solução exata, termo fonte e adequa o tensor de
        % permeabilidade
        for ielem=1:size(elem,1)
            x=centelem(ielem,1);
            y=centelem(ielem,2);
            if centelem(ielem,1)<0.5
                u(ielem,1)=c1*x^2+d1*y^2;
                kmap(ielem,1:5)=[ielem 50 0 0 1];
                elem(ielem,5)=ielem;
                fonte(ielem,1)=-100*c1-2*d1;
                normKmap(ielem,1)=norm([50 0; 0 1]);
            else
                u(ielem,1)=1+b2*x+c2*x^2+d2*y^2;
                kmap(ielem,1:5)=[ielem 1 0 0 10];
                elem(ielem,5)=ielem;
                fonte(ielem,1)=-2*c2-20*d2;
                normKmap(ielem,1)=norm([1 0; 0 10]);
            end
        end
        % calculo das velocidades
        R=[0 -1 0; 1 0 0; 0 0 0];
        for i=1:size(bedge,1)
            
            IJ=coord(bedge(i,2),:)-coord(bedge(i,1),:);
            norma=norm(IJ);
            nij=R*IJ'/norma;
            auxpoint=(coord(bedge(i,2),:)+coord(bedge(i,1),:))*0.5;
            x=auxpoint(1,1);
            y=auxpoint(1,2);
            if x<0.5
                a= [-100*c1*x -2*d1*y 0];
            else
                a= [-b2-2*c2*x -20*d2*y 0];
            end
            F(i,1)=dot(a,nij');
            
        end
        
        for i=1:size(inedge,1)
            
            IJ=coord(inedge(i,2),:)-coord(inedge(i,1),:);
            norma=norm(IJ);
            nij=R*IJ'/norma;
            auxpoint=(coord(inedge(i,2),:)+coord(inedge(i,1),:))*0.5;
            x=auxpoint(1,1);
            y=auxpoint(1,2);
            if x<0.5
                a= [-100*c1*x -2*d1*y 0];
            else
                a= [-b2-2*c2*x -20*d2*y 0];
            end
            F(i+size(bedge,1),1)=dot(a,nij');
            
        end
        vel=F;
        K=kmap;
    case 'FriisEdwards'
        u=0;
        fonte=zeros(size(elem,1),1);
        %[element1]=elementsneigbormiddleface;
        % 253   five_spotMALHAV1
        % 5178  five_spotMALHAV2
        % 21517 five_spotMALHAV3
        % 20064 five_spotMALHAV4
        % 9985 MESH_06TRIANGULARDISTORTED_96PROBLEM431
        
        fonte(9985,1)=248000000;
        
        k = [3000 0; 0 1];
        %Fill "R"
        theta=-pi/6;
        R(1,1) = cos(theta);
        R(1,2) = sin(theta);
        R(2,1) = -R(1,2);
        R(2,2) = R(1,1);
        %Fill "k" turning the tensor
        A=inv(R);
        k = A*k*R;
        %Buld "kmap" again
        kmap = [1 k(1,1) k(1,2) k(2,1) k(2,2)];
        elem(:,5)=1;
        K=kmap;
    case 'shenyuan16'
        
        for i=1:size(elem,1)
            x= centelem(i,1);
            y=centelem(i,2);
            k = [1+2*x^2+y^2 0; 0 1+x^2+2*y^2];
            %Fill "R"
            theta=5*pi/12;
            R(1,1) = cos(theta);
            R(1,2) = sin(theta);
            R(2,1) = -R(1,2);
            R(2,2) = R(1,1);
            %Fill "k" turning the tensor
            A=inv(R);
            k = A*k*R;
            %Buld "kmap" again
            kmap(i,1:5) = [i k(1,1) k(1,2) k(2,1) k(2,2)];
            % valido
            t1=2.46711*(x^2-y^2)*cos(pi*x)*cos(pi*y)-pi^2*(1+1.0669*x^2+1.93289*y^2)*sin(pi*x)*sin(pi*y)+...
                1.57061*x*sin(pi*x)*cos(pi*y)+6.70353*x*cos(pi*x)*sin(pi*y);
            
            t2=2.46711*(x^2-y^2)*cos(pi*x)*cos(pi*y)-pi^2*(1+1.93289*x^2+1.0669*y^2)*sin(pi*x)*sin(pi*y)+...
                6.70353*y*sin(pi*x)*cos(pi*y)-1.57061*y*cos(pi*x)*sin(pi*y);
            fonte(i,1)=-(t1+t2)*elemarea(i,1) ;
            
            u(i,1)=sin(pi*x)*sin(pi*y);
            elem(i,5)=i;
        end
        for iface=1:size(bedge,1)+size(inedge,1)
            R1=[0 -1 0; 1 0 0; 0 0 0];
            if iface< size(bedge,1) || iface==size(bedge,1)
                v1=bedge(iface,1);
                v2=bedge(iface,2);
                IJ=coord(v2,:)-coord(v1,:);
                norma=norm(IJ);
                nij=R1*IJ'/norma;
                auxpoint=(coord(v2,:)+coord(v1,:))*0.5;
                x=auxpoint(1,1);
                y=auxpoint(1,2);
                t0= pi*(1+1.0670*x^2+1.9330*y^2)*cos(pi*x)*sin(pi*y)+pi*(x^2-y^2)*0.25*sin(pi*x)*cos(pi*y);
                t01=pi*(1+1.9330*x^2+1.0670*y^2)*sin(pi*x)*cos(pi*y)+pi*(x^2-y^2)*0.25*cos(pi*x)*sin(pi*y);
                
                % anterior
                
                %t0= pi*(1+1.0669*x^2+1.93289*y^2)*cos(pi*x)*sin(pi*y)+pi*(x^2-y^2)*0.24997*sin(pi*x)*cos(pi*y);
                %t01=pi*(1+1.93289*x^2+1.0669*y^2)*sin(pi*x)*cos(pi*y)+pi*(x^2-y^2)*0.24997*cos(pi*x)*sin(pi*y);
                
                
                a=[-t0, -t01, 0];
            else
                v1=inedge(iface-size(bedge,1),1);
                v2=inedge(iface-size(bedge,1),2);
                IJ=coord(v2,:)-coord(v1,:);
                norma=norm(IJ);
                nij=R1*IJ'/norma;
                auxpoint=(coord(v2,:)+coord(v1,:))*0.5;
                x=auxpoint(1,1);
                y=auxpoint(1,2);
                t0= pi*(1+1.0670*x^2+1.9330*y^2)*cos(pi*x)*sin(pi*y)+pi*(x^2-y^2)*0.25*sin(pi*x)*cos(pi*y);
                t01=pi*(1+1.9330*x^2+1.0670*y^2)*sin(pi*x)*cos(pi*y)+pi*(x^2-y^2)*0.25*cos(pi*x)*sin(pi*y);
                
                % anterior
                %t0= pi*(1+1.0669*x^2+1.93289*y^2)*cos(pi*x)*sin(pi*y)+pi*(x^2-y^2)*0.24997*sin(pi*x)*cos(pi*y);
                %t01=pi*(1+1.93289*x^2+1.0669*y^2)*sin(pi*x)*cos(pi*y)+pi*(x^2-y^2)*0.24997*cos(pi*x)*sin(pi*y);
                
                
                a=[-t0, -t01, 0];
                
            end
            F(iface,1)=dot(a,nij');
            
        end
        vel=F;
        K=kmap;
    case 'gaowu5'
        fonte=0;
        elem(:,5)=1;
        u=0;
        %Initialize "R":
        R = zeros(2);
        k = [1 0; 0 1e-4];
        %Fill "R"
        R(1,1) = cosd(40);
        R(1,2) = sind(40);
        R(2,1) = -R(1,2);
        R(2,2) = R(1,1);
        %Fill "k" turning the tensor
        A=inv(R);
        k = A*k*R;
        %Buld "kmap" again
        kmap = [1 k(1,1) k(1,2) k(2,1) k(2,2)];
        K=kmap;
        
    case 'edqueiroz'
        % problema artigo quiroz et al.
        epsilon=1e3;
        for ielem=1:size(elem,1)
            if centelem(ielem,1)<0.5 ||  centelem(ielem,1)==0.5
                theta=0.5*pi;
                k=[100 0;0 0.01];
                karot=[ cos(theta) -sin(theta); sin(theta) cos(theta)]*k*...
                    [cos(theta) sin(theta); -sin(theta) cos(theta)];
                kmap(ielem,1:5)=[ielem karot(1,1) karot(1,2) karot(2,1) karot(2,2)];
                elem(ielem,5)=ielem;
            else
                x1= centelem(ielem,1)+1e-3;
                y1= centelem(ielem,2)+1e-3;
                
                kmap(ielem,1:5)=[ielem y1^2+epsilon*x1^2 -(1-epsilon)*x1*y1 -(1-epsilon)*x1*y1 x1^2+epsilon*y1^2 ];
                elem(ielem,5)=ielem;
                
            end
        end
        K=kmap;
        fonte=0;
        u=0;
        %% use se deseja resolver o problema de Queiroz et al 2014
        % unstructured mesh com furo reordenando o sentido da fronterira no
        % contorno interior
        % unstructured mesh com furo reordenando o sentido da fronterira no
        % contorno interior
        %          x=bedge(71:78,1); % UTILIZE Benchmark23_3_18_18.msh
        %          y=bedge(71:78,2);
        %          bedge(71:78,1)=y;
        %          bedge(71:78,2)=x;
        %          bedge(71:78,4:5)=102; % 18x18
        %          bcflag(2,1)=102;
        %          bcflag(2,2)=2;
        %=====================================
        %          x=bedge(73:80,1); % UTILIZE MeshTri18sym.msh
        %          y=bedge(73:80,2);
        %          bedge(73:80,1)=y;
        %          bedge(73:80,2)=x;
        %          bedge(73:80,4:5)=102; % 18x18
        %          bcflag(2,1)=102;
        %          bcflag(2,2)=2;
        
        %=====================================
        %
        x=bedge(135:150,1); % UTILIZE Benchmark23_3_FINA36.msh
        y=bedge(135:150,2);
        bedge(135:150,1)=y;
        bedge(135:150,2)=x;
        bedge(135:150,4:5)=102; % 36x36
        bcflag(2,1)=102;
        bcflag(2,2)=2;
        %===============================================================
        %  x=bedge(145:160,1); % UTILIZE MeshTri36sym.msh
        %  y=bedge(145:160,2);
        %  bedge(145:160,1)=y;
        %  bedge(145:160,2)=x;
        %  bedge(145:160,4:5)=102; % 36x36
        %  bcflag(2,1)=102;
        %  bcflag(2,2)=2;
        
        %===============================================================
        % unstructured mesh com furo reordenando o sentido da fronterira no
        % contorno interior
        %         x=bedge(289:320,1);
        %         y=bedge(289:320,2);
        %         bedge(289:320,1)=y;
        %         bedge(289:320,2)=x;
        %         bedge(289:320,4:5)=102; % benchmark23_3 72x72
        %         bcflag(2,1)=102;
        %         bcflag(2,2)=2;
        
        % unstructured mesh com furo reordenando o sentido da fronterira no
        % contorno interior
        %  x=bedge(577:640,1);
        %  y=bedge(577:640,2);
        %  bedge(577:640,1)=y;
        %  bedge(577:640,2)=x;
        %  bedge(577:640,4:5)=102; % benchmark23_3 144x144
        %  bcflag(2,1)=102;
        %  bcflag(2,2)=2;
        
    case 'homogeneo'
        %% meio isotropico homogeneo
        for i=1:size(centelem,1)
            
            kmap(1,1:5)=[1 1 0 0 1];
            elem(i,5)=1;
            
            u(i,1)=1-centelem(i,1);
            fonte=0;
        end
        
        for iface=1:size(bedge,1)+size(inedge,1)
            R1=[0 -1 0; 1 0 0; 0 0 0];
            if iface< size(bedge,1) || iface==size(bedge,1)
                v1=bedge(iface,1);
                v2=bedge(iface,2);
                IJ=coord(v2,:)-coord(v1,:);
                norma=norm(IJ);
                nij=R1*IJ'/norma;
                
                a=[1, 0, 0];
            else
                v1=inedge(iface-size(bedge,1),1);
                v2=inedge(iface-size(bedge,1),2);
                IJ=coord(v2,:)-coord(v1,:);
                norma=norm(IJ);
                nij=R1*IJ'/norma;
                
                a=[1, 0, 0];
                
            end
            F(iface,1)=dot(a,nij');
            
        end
        vel=F;
        K=kmap;
    case 'heterogeneo'
        for i=1:size(elem,1)
            if centelem(i,1)<0.5
                kmap(1,1:5)=[1 1 0 0 1];
                elem(i,5)=1;
                u(i,1)=(4/3)*centelem(i,1);
            else
                
                kmap(2,1:5)=[1 2 0 0 2];
                
                elem(i,5)=2;
                u(i,1)=(2/3)*centelem(i,1)+(1/3);
            end
            fonte=0;
        end
        for iface=1:size(bedge,1)+size(inedge,1)
            R1=[0 -1 0; 1 0 0; 0 0 0];
            if iface< size(bedge,1) || iface==size(bedge,1)
                v1=bedge(iface,1);
                v2=bedge(iface,2);
                IJ=coord(v2,:)-coord(v1,:);
                norma=norm(IJ);
                nij=R1*IJ'/norma;
                if IJ(1,1)<0.5
                    
                    a=[-4/3, 0, 0];
                else
                    a=[-4/3, 0, 0];
                end
            else
                v1=inedge(iface-size(bedge,1),1);
                v2=inedge(iface-size(bedge,1),2);
                IJ=coord(v2,:)-coord(v1,:);
                norma=norm(IJ);
                nij=R1*IJ'/norma;
                
                if IJ(1,1)<0.5
                    
                    a=[-4/3, 0, 0];
                else
                    a=[-4/3, 0, 0];
                end
                
            end
            F(iface,1)=dot(a,nij');
            
        end
        vel=F;
        K=kmap;
    case 'crumpton'
        
        alfa=1000;
        u=zeros(size(elem,1),1);
        for ielem=1:size(elem,1)
            if centelem(ielem,1)<0 || centelem(ielem,1)==0
                u(ielem)=(2*sin(centelem(ielem,2))+cos(centelem(ielem,2)))*alfa*centelem(ielem,1)+sin(centelem(ielem,2))+3*alfa;
                kmap(1,1:5)=[1 1 0 0 1];
                elem(ielem,5)=1;
                fonte(ielem,1)=(2*sin(centelem(ielem,2))+cos(centelem(ielem,2)))*alfa*centelem(ielem,1)+sin(centelem(ielem,2));
            else
                u(ielem)=exp(centelem(ielem,1))*sin(centelem(ielem,2))+3*alfa;
                kmap(2,1:5)=[2 alfa*2 alfa*1 alfa*1 alfa*2];
                elem(ielem,5)=2;
                fonte(ielem,1)=-2*alfa*exp(centelem(ielem,1))*cos(centelem(ielem,2));
            end
        end
        
        % calculo dos fluxos analiticos
        R=[0 -1 0; 1 0 0; 0 0 0];
        for i=1:size(bedge,1)
            
            IJ=coord(bedge(i,2),:)-coord(bedge(i,1),:);
            norma=norm(IJ);
            nij=R*IJ'/norma;
            auxpoint=(coord(bedge(i,2),:)+coord(bedge(i,1),:))*0.5;
            x=auxpoint(1,1);
            y=auxpoint(1,2);
            if x<0 || x==0
                a= [-(2*sin(y)+ cos(y))*alfa -(2*cos(y)- sin(y))*alfa*x-cos(y) 0];
            else
                a= alfa*exp(x)*[-(2*sin(y)+ cos(y)) -(sin(y)+ 2*cos(y)) 0];
            end
            F(i,1)=dot(a,nij');
            
        end
        
        for i=1:size(inedge,1)
            
            IJ=coord(inedge(i,2),:)-coord(inedge(i,1),:);
            norma=norm(IJ);
            nij=R*IJ'/norma;
            auxpoint=(coord(inedge(i,2),:)+coord(inedge(i,1),:))*0.5;
            x=auxpoint(1,1);
            y=auxpoint(1,2);
            if x<0 || x==0
                a= [-(2*sin(y)+ cos(y))*alfa -(2*cos(y)- sin(y))*alfa*x-cos(y) 0];
            else
                a= alfa*exp(x)*[-(2*sin(y)+ cos(y)) -(sin(y)+ 2*cos(y)) 0];
            end
            F(i+size(bedge,1),1)=dot(a,nij');
            
        end
        K=kmap;
        vel=F;
    case 'gaowu1'
        for i=1:size(elem,1)
            x=centelem(i,1);
            y=centelem(i,2);
            if x<0.5 || x==0.5
                fonte(i,1)=-exp(-62.8319*((-0.5+y)^2))*(127304-782890*y+(1.55068*(10^6))*...
                    (y^2)-992201*(y^3)+x*(-222452+(1.43303*10^6)*y-(2.96871*10^6)*(y^2)+...
                    (1.9844*10^6)*(y^3)));
                % obtido no wolfram
                
                %                fonte(i,1)=-exp(-62.8319*((y-0.5)^2))*(128821-789183*y+(1.557*(10^6))*...
                %                     (y^2)-992204*(y^3)+x*(-222452+(1.43303*10^6)*y-(2.96871*10^6)*(y^2)+...
                %                     (1.9844*10^6)*(y^3)));
                
                u(i,1)=(1+(x-0.5)*(0.1+8*pi*(y-0.5)))*exp(-20*pi*((y-0.5)^2));
                elem(i,5)=1;
                kmap(1,:)=[1 10 2 2 5];
            else
                fonte(i,1)=(-exp(x-62.8319*((y-0.5)^2))*(2318.87-9577.95*y+9577.95*(y^2)));
                % obtido por wolfram
                %fonte(i,1)=-exp(x-62.8319*((y-0.5)^2))*(23.789231-95.7796*y+95.7796*(y^2));
                u(i,1)=exp(x-0.5)*exp(-20*pi*((y-0.5)^2));
                
                elem(i,5)=2;
                kmap(2,:)=[2 1 0 0 1];
            end
        end
        K=kmap;
    case 'gaowu2'
        c=centelem;
        R=[0 -1 0; 1 0 0; 0 0 0];
        for i=1:size(elem,1)
            x=c(i,1);
            y=c(i,2);
            % valido
            fonte(i,1)= (-0.25*(-3*(-1 + c(i,1))^2*(-1 + c(i,2)) + cos((-1 + c(i,1))*(-1 + c(i,2)))* csc(1)) -  0.25*(-1 + c(i,2))*(-3*(-1 + c(i,1))^2 - (-1 + c(i,1))*csc(1)*sin((-1 + c(i,1))*(-1 + c(i,2)))) - ...
                0.75*(-2*(-1 + c(i,1))^3 - (-1 + c(i,1))^2*csc(1)*sin((-1 + c(i,1))*(-1 + c(i,2)))) - 0.75*(-1 + c(i,2))*(-6*(-1 + c(i,1))*(-1 + c(i,2)) - (-1 + c(i,2))*csc(1)* sin((-1 + c(i,1))*(-1 + c(i,2)))) - ...
                0.25*(-6*(-1 + c(i,1))^2 *(-1 + c(i,2)) +cos((-1 + c(i,1))*(-1 + c(i,2)))*csc(1) - (-1 + c(i,1))*(-1 + c(i,2))*csc(1)*sin((-1 + c(i,1))*(-1 + c(i,2)))));
            %             fonte(i,1)=-0.25*(-3*((x - 1)^2)*(y - 1) + ...
            %                     cos((x - 1)*(y - 1))*csc(1)) - ...
            %                     0.25*(y - 1)*(-3*((x - 1)^2) - ...
            %                     (x - 1)*csc(1)*sin((x - 1)*(y - 1))) - 0.75*(-2*((x - 1)^3) - ...
            %                     ((x - 1)^2)*csc(1)*sin((x - 1)*(y - 1))) - ...
            %                     0.75*(y - 1)*(-6*(x - 1)*(y - 1) - ...
            %                     (y - 1)*csc(1)*sin((x - 1)*(y - 1))) - ...
            %                     0.25*(-6*((x - 1)^2)*(y - 1) + cos((x - 1)*(y - 1))*csc(1) - ...
            %                     (x - 1)*(y - 1)*csc(1)*sin((x - 1)*(y - 1)));
            
            x= centelem(i,1);
            y=centelem(i,2);
            u(i,1)=0.5*((sin((1-x)*(1-y))/(sin(1)))+(1-x)^3*(1-y)^2);
            elem(i,5)=1;
        end
        
        for iface=1:size(bedge,1)+size(inedge,1)
            
            if iface< size(bedge,1) || iface==size(bedge,1)
                v1=bedge(iface,1);
                v2=bedge(iface,2);
                IJ=coord(v2,:)-coord(v1,:);
                norma=norm(IJ);
                nij=R*IJ'/norma;
                auxpoint=(coord(v2,:)+coord(v1,:))*0.5;
                x=auxpoint(1,1);
                y=auxpoint(1,2);
                
                dpdx=((1-y)/2)*(((-cos((1-x)*(1-y)))/sin(1))-3*(1-x)^2*(1-y));
                dpdy=((1-x)/2)*(((-cos((1-x)*(1-y)))/sin(1))-2*(1-x)^2*(1-y));
                a=-[1.5*dpdx+0.5*dpdy 0.5*dpdx+1.5*dpdy 0];
                
            else
                v1=inedge(iface-size(bedge,1),1);
                v2=inedge(iface-size(bedge,1),2);
                IJ=coord(v2,:)-coord(v1,:);
                norma=norm(IJ);
                nij=R*IJ'/norma;
                auxpoint=(coord(v2,:)+coord(v1,:))*0.5;
                x=auxpoint(1,1);
                y=auxpoint(1,2);
                
                
                dpdx=((1-y)/2)*(((-cos((1-x)*(1-y)))/sin(1))-3*(1-x)^2*(1-y));
                dpdy=((1-x)/2)*(((-cos((1-x)*(1-y)))/sin(1))-2*(1-x)^2*(1-y));
                a=-[1.5*dpdx+0.5*dpdy 0.5*dpdx+1.5*dpdy 0];
                
            end
            F(iface,1)=dot(a,nij');
            
        end
        
        kmap=[1 1.5 0.5 0.5 1.5];
        K=kmap;
        vel=F;
    case 'crumptonhyman'
        
        for i=1:size(elem,1)
            x= centelem(i,1);
            y=centelem(i,2);
            % valido
            fonte(i,1)= -2*(1+x^2+x*y+y^2)*exp(x*y);
            u(i,1)=exp(x*y);
            elem(i,5)=1;
        end
        kmap=[1 2 1 1 2];
    case 'gaowu3'
        alfa=1000;
        for i=1:size(elem,1)
            x=centelem(i,1);
            y=centelem(i,2);
            r=exp(-20*pi*((x-0.5)^2+(y-0.5)^2));
            %% a mano
            %             fonte(i,1)=-r*(-40*pi*(3*alfa-1)*(x*(x-0.5)+y*(y-0.5))+...
            %                 3200*pi^2*x*y*(x-0.5)*(y-0.5)*(alfa-1)+...
            %                 (alfa*x^2+y^2)*(-40*pi+1600*pi^2*(x-0.5)^2)+...
            %                 (x^2+alfa*y^2)*(-40*pi+1600*pi^2*(y-0.5)^2));
            %% wolfram alpha
            
            parte1=(y^2+alfa*x^2)*r*(- 125.664 + 15791.4*(x-0.5)^2)+...
                r*15791.4*(alfa-1)*(y-0.5)*y*x*(x-0.5)- 251.328*alfa*x*(x-0.5)*r-...
                125.664*(alfa-1)*(y-0.5)*y*r ;
            
            parte2=r*(x^2+alfa*y^2)*(-125.664 + 15791.4*(y-0.5)^2)+...
                r*15791.4*(alfa-1)*(x-0.5)*x*y*(y-0.5) - 251.328*alfa*y*(y-0.5)*r-...
                125.664*(alfa-1)*(x-0.5)*x*r;
            
            fonte(i,1)=-parte1-parte2;
            u(i,1)=exp(-20*pi*((x-0.5)^2 + (y-0.5)^2))+1e-3;
            elem(i,5)=i;
            kmap(i,1:5)=[i alfa*x^2+y^2 (alfa-1)*x*y (alfa-1)*x*y x^2+alfa*y^2];
        end
        
    case 'lepotier'
        R=[0 1 0; -1 0 0; 0 0 0];
        alfa=1e-2;
        for i=1:size(elem,1)
            x=centelem(i,1);
            y=centelem(i,2);
            uu=sin(pi*x)*sin(pi*y);
            r=cos(pi*x)*sin(pi*y);
            s=sin(pi*x)*cos(pi*y);
            x1=x+alfa;
            y1=y+alfa;
            %% no artigo de Le potier
            fonte(i,1)=-elemarea(i,1)*(-uu*((1-alfa)*pi^2*(x1^2+y1^2))-r*((1-3*alfa)*pi*x1)-s*((1-3*alfa)*pi*y1)-...
                cos(pi*x)*cos(pi*y)*(2*pi^2*(1-alfa)*x1*y1));
            k=1;
            u(i,1)=sin(pi*x)*sin(pi*y);
            % in the article Zhang Kobaise 2019, consider:
            %             fonte(i,1)=elemarea(i,1)*(pi/(x^2+y^2))*(-(alfa - 1)*x*...
            %                 cos(pi*x)*(2*pi*y*...
            %                 cos(pi*y) + sin(pi*y)) + ...
            %                 sin(pi*x)*(-(alfa - 1)*y*...
            %                 cos(pi*y) + (1 + alfa)*pi*((x^2) + (y^2))*...
            %                 sin(pi*y)));
            %             k = (1/((x^2) + (y^2)));
            %             u(i,1)=sin(pi*x)*sin(pi*y)+1;
            elem(i,5)=i;
            kmap(i,1:5)=[i k*(alfa*x^2+y^2) k*((alfa-1)*x*y) k*((alfa-1)*x*y) k*(x^2+alfa*y^2)];
        end
        K=kmap;
        for iface=1:size(bedge,1)+size(inedge,1)
            
            if iface< size(bedge,1) || iface==size(bedge,1)
                v1=bedge(iface,1);
                v2=bedge(iface,2);
                IJ=coord(v2,:)-coord(v1,:);
                norma=norm(IJ);
                nij=R*IJ'/norma;
                auxpoint=(coord(v2,:)+coord(v1,:))*0.5;
                x=auxpoint(1,1);
                y=auxpoint(1,2);
                
            else
                v1=inedge(iface-size(bedge,1),1);
                v2=inedge(iface-size(bedge,1),2);
                IJ=coord(v2,:)-coord(v1,:);
                norma=norm(IJ);
                nij=R*IJ'/norma;
                auxpoint=(coord(v2,:)+coord(v1,:))*0.5;
                x=auxpoint(1,1);
                y=auxpoint(1,2);
            end
            
            dpdx=pi*cos(pi*x)*sin(pi*y);
            dpdy=pi*sin(pi*x)*cos(pi*y);
            k11=(y^2+alfa*x^2)/(x^2 + y^2);
            k12=((alfa-1)*x*y)/(x^2 + y^2);
            k22=(alfa*y^2+x^2)/(x^2 + y^2);
            
            a=-[k11*dpdx+k12*dpdy  k12*dpdx+k22*dpdy 0];
            
            F(iface,1)=dot(a,(R*IJ')');
            
        end
        vel=F;
    case 'gaowu4'
        alfa=1e-6;
        for i=1:size(elem,1)
            x=centelem(i,1);
            y=centelem(i,2);
            uu=sin(pi*x)*sin(pi*y);
            r=cos(pi*x)*sin(pi*y);
            s=sin(pi*x)*cos(pi*y);
            d=x^2+y^2;
            %% obtido por wolfram alpha1
            %             part1= (1/d)*(-pi^2*uu*(alfa*x^2+y^2)+ pi^2*(alfa-1)*x*y*cos(pi*x)*cos(pi*y)+ ...
            %                 pi*(alfa-1)*y*s+ 2*pi*alfa*x*r)- (1/d^2)*(2*x*(pi*r*(alfa*x^2+y^2)+pi*(alfa-1)*x*y*s));
            %
            %             part2=(1/d)*(-pi^2*uu*(alfa*y^2+x^2)+ pi^2*(alfa-1)*x*y*cos(pi*x)*cos(pi*y)+ ...
            %                 2*pi*alfa*y*s+ pi*(alfa-1)*x*r)-(1/d^2)*(2*y*(pi*s*(alfa*y^2+x^2)+pi*(alfa-1)*x*y*r));
            %% wolfram alpha2
            
            part1= (1/d)*(-pi^2*uu*(alfa*x^2+y^2)- pi^2*(alfa-1)*x*y*uu+ ...
                pi*(alfa-1)*x*s+ 2*pi*alfa*x*r)- (1/d^2)*(2*x*pi*(r*(alfa*x^2+y^2)+(alfa-1)*y^2*s));
            
            
            
            part2=(1/d)*(-pi^2*uu*(alfa*y^2+x^2)- pi^2*(alfa-1)*x*y*uu+ ...
                2*pi*alfa*y*s+ pi*(alfa-1)*y*r) -   (1/d^2)*(2*y*pi*(s*(alfa*y^2+x^2)+(alfa-1)*x^2*r));
            
            fonte(i,1)=-part1-part2;
            
            u(i,1)=sin(pi*x)*sin(pi*y);
            elem(i,5)=i;
            kmap(i,1:5)=[i (alfa*x^2+y^2)/d ((alfa-1)*x*y)/d ((alfa-1)*x*y)/d (x^2+alfa*y^2)/d];
            
        end
        K=kmap;
    case 'lipnikov1'
        
        for i=1:size(elem,1)
            x=centelem(i,1);
            y=centelem(i,2);
            if x<0.5 || x==0.5
                fonte(i,1)=4;
                u(i,1)=1-2*y^2+4*x*y+6*x+2*y;
                elem(i,5)=1;
                kmap(1,1:5)=[1 1 0 0 1];
            else
                fonte(i,1)=-5.6 ;
                u(i,1)=-2*y^2+1.6*x*y-0.6*x+3.2*y+4.3;
                elem(i,5)=2;
                kmap(2,1:5)=[2 10 3 3 1];
            end
            
        end
        
    case 'lipnikov2'
        u=0;
        fonte=0;
        k = [1000 0; 0 1];
        %Fill "R"
        %theta=5*pi/6;
        theta=67.5;
        R(1,1) = cosd(theta);
        R(1,2) = sind(theta);
        R(2,1) = -R(1,2);
        R(2,2) = R(1,1);
        %Fill "k" turning the tensor
        A=inv(R);
        k = A*k*R;
        %Buld "kmap" again
        kmap = [1 k(1,1) k(1,2) k(2,1) k(2,2)];
        elem(:,5)=1;
        K=kmap;
        
    case 'lipnikov3'
        u=0;
        fonte=0;
        k = [100000 0; 0 1];
        %Fill "R"
        %theta=pi/4;
        theta=-pi/6;
        R(1,1) = cos(theta);
        R(1,2) = sin(theta);
        R(2,1) = -R(1,2);
        R(2,2) = R(1,1);
        %Fill "k" turning the tensor
        A=inv(R);
        k = A*k*R;
        %Buld "kmap" again
        kmap = [1 k(1,1) k(1,2) k(2,1) k(2,2)];
        elem(:,5)=1;
        K=kmap;
    case 'herbinhubert'
        %Initialize a parameter
        epsilon = 1e3;
        theta = 0.5*pi;
        k = [100 0; 0 0.01];
        %Rotate the tensor in "theta"
        Krot = [cos(theta) -sin(theta); sin(theta) cos(theta)]*k*...
            [cos(theta) sin(theta); -sin(theta) cos(theta)];
        
        %Initialize "kmap"
        kmap = zeros(size(centelem,1),5);
        for i = 1:size(centelem,1)
            %Define "x" and "y"
            x = centelem(i,1);
            y = centelem(i,2);
            %Choose the tensor according "x" position
            if x < 0.5
                kmap (i,:) = [i Krot(1,1) Krot(1,2) Krot(2,1) Krot(2,2)];
            else
                %Define "x1" and "y1"
                x1 = x + 1e-3;
                y1 = y + 1e-3;
                %Definition of permeability components
                k(1,1) = ((y1^2) + epsilon*(x1^2));
                k(1,2) = -(1 - epsilon)*(x1*y1);
                k(2,1) = -(1 - epsilon)*(x1*y1);
                k(2,2) = ((x1^2) + epsilon*(y1^2));
                %Build "kmap"
                kmap(i,:) = [i k(1,1) k(1,2) k(2,1) k(2,2)];
            end  %End of IF
        end  %End of FOR
    case 'guangwei'
        u=0;
        alfa=5*1e-3;
        for i=1:size(centelem,1)
            x=centelem(i,1);
            y=centelem(i,2);
            kmap(i,1:5)=[i (y^2+alfa*x^2)+alfa -(1-alfa)*x*y -(1-alfa)*x*y  y^2*alfa+x^2+alfa];
            elem(i,5)=i;
            if (x>3/8 || x==3/8)&&(x<5/8 || x==5/8)&&(y>3/8 || y==3/8)&&(y<5/8 || y==5/8)
                fonte(i,1)=1;
            else
                fonte(i,1)=0;
            end
        end
        K=kmap;
    case 'guangwei1'
        
        u=0;
        %Initialize "R":
        R = zeros(2);
        k = [10000 0; 0 1];
        %Fill "R"
        R(1,1) = cos(pi/6);
        R(1,2) = sin(pi/6);
        R(2,1) = -R(1,2);
        R(2,2) = R(1,1);
        %Fill "k" turning the tensor
        A=inv(R);
        k = A*k*R;
        %Buld "kmap" again
        kmap = [1 k(1,1) k(1,2) k(2,1) k(2,2)];
        for i=1:size(centelem,1)
            x=centelem(i,1);
            y=centelem(i,2);
            elem(i,5)=1;
            %             if (7/18<centelem(i,1) || 7/18==centelem(i,1))&& ...
            %                     (centelem(i,1)<11/18 || 11/18==centelem(i,1))...
            %                     && (7/18<centelem(i,2) || centelem(i,2)==7/18) && ...
            %                     (centelem(i,2)<11/18 || centelem(i,2)==11/18)
            if (50/101<x || 50/101==x)&& (x<51/101 || 51/101==x)...
                    && (50/101<y || y==50/101) && (y<51/101 || y==51/101)
                fonte(i,1)=1021;
            else
                fonte(i,1)=0;
            end
        end
        K=kmap;
    case 'gaowu7'
        auxbedge1=bedge(:,1);
        auxbedge2=bedge(:,2);
        
        bedge(:,1)=auxbedge2;
        bedge(:,2)=auxbedge1;
        
        u=0;
        alfa=0.2;
        for i=1:size(centelem,1)
            x=centelem(i,1);
            y=centelem(i,2);
            phi1=y-alfa*(x-0.5)-0.475;
            phi2=phi1-0.05;
            theta=atand(alfa);
            if phi1<0 || phi1==0
                elem(i,5)=i;
                %Initialize "R":
                R = zeros(2);
                k = [1 0; 0 0.1];
                %Fill "R"
                R(1,1) = cosd(theta);
                R(1,2) = sind(theta);
                R(2,1) = -R(1,2);
                R(2,2) = R(1,1);
                %Fill "k" turning the tensor
                A=inv(R);
                k = A*k*R;
                %Buld "kmap" again
                kmap(i,1:5) = [i k(1,1) k(1,2) k(2,1) k(2,2)];
                fonte(i,1)=0;
            elseif phi1>0 && phi2<0
                elem(i,5)=i;
                %Initialize "R":
                R = zeros(2);
                k = [100 0; 0 10];
                %Fill "R"
                R(1,1) = cosd(theta);
                R(1,2) = sind(theta);
                R(2,1) = -R(1,2);
                R(2,2) = R(1,1);
                %Fill "k" turning the tensor
                A=inv(R);
                k = A*k*R;
                %Buld "kmap" again
                kmap(i,1:5) = [i k(1,1) k(1,2) k(2,1) k(2,2)];
                fonte(i,1)=0;
            elseif phi2>0 || phi2==0
                elem(i,5)=i;
                %Initialize "R":
                R = zeros(2);
                k = [1 0; 0 0.1];
                %Fill "R"
                R(1,1) = cosd(theta);
                R(1,2) = sind(theta);
                R(2,1) = -R(1,2);
                R(2,2) = R(1,1);
                %Fill "k" turning the tensor
                A=inv(R);
                k = A*k*R;
                %Buld "kmap" again
                kmap(i,1:5) = [i k(1,1) k(1,2) k(2,1) k(2,2)];
                fonte(i,1)= 0;
            end
            
            u(i,1)=-x-alfa*y;
            KK=kmap(i,2:5);
            r=[KK(1,1) KK(1,2);
                KK(1,3) KK(1,4)];
            normKmap(i,1)=norm(r);
        end
        K=kmap;
        
    case 'gaowu8'
        auxbedge1=bedge(:,1);
        auxbedge2=bedge(:,2);
        
        bedge(:,1)=auxbedge2;
        bedge(:,2)=auxbedge1;
        
        u=0;
        alfa=0.2;
        fonte=0;
        for i=1:size(centelem,1)
            x=centelem(i,1);
            y=centelem(i,2);
            phi1=y-alfa*(x-0.5)-0.475;
            phi2=phi1-0.05;
            
            if phi1<0
                elem(i,5)=i;
                %Buld "kmap" again
                kmap(i,1:5) = [i 1 0 0 1];
                u(i,1)=-phi1;
            elseif phi1>0 && phi2<0
                elem(i,5)=i;
                %Buld "kmap" again
                kmap(i,1:5) = [i 0.01 0 0 0.01];
                u(i,1)=-100*phi1;
                
            elseif phi2>0
                elem(i,5)=i;
                
                %Buld "kmap" again
                kmap(i,1:5) = [i 1 0 0 1];
                u(i,1)=-phi2-5;
            end
            KK=kmap(i,2:5);
            r=[KK(1,1) KK(1,2);
                KK(1,3) KK(1,4)];
            normKmap(i,1)=norm(r);
        end
        for iface=1:size(bedge,1)+size(inedge,1)
            R=[0 -1 0; 1 0 0; 0 0 0];
            if iface< size(bedge,1) || iface==size(bedge,1)
                v1=bedge(iface,1);
                v2=bedge(iface,2);
                IJ=coord(v2,:)-coord(v1,:);
                norma=norm(IJ);
                nij=R*IJ'/norma;
                auxpoint=(coord(v2,:)+coord(v1,:))*0.5;
                x=auxpoint(1,1);
                y=auxpoint(1,2);
                phi1=y-alfa*(x-0.5)-0.475;
                phi2=phi1-0.05;
                if phi1<0
                    
                    %Initialize "R":
                    
                    k = [1 0; 0 1];
                    
                    %Buld "kmap" again
                    KK = [k(1,1) k(1,2) 0;k(2,1) k(2,2) 0; 0 0 0];
                    a=-KK*[alfa;-1;0];
                elseif phi1>0 && phi2<0
                    
                    %Initialize "R":
                    
                    k = [0.01 0; 0 0.01];
                    
                    %Buld "kmap" again
                    KK= [k(1,1) k(1,2) 0; k(2,1) k(2,2) 0; 0 0 0];
                    a=-KK*[100*alfa;-100;0];
                elseif phi2>0
                    
                    %Initialize "R":
                    
                    k = [1 0; 0 1];
                    
                    %Buld "kmap" again
                    KK = [k(1,1) k(1,2) 0; k(2,1) k(2,2) 0; 0 0 0];
                    a=-KK*[alfa;-1;0];
                end
            else
                v1=inedge(iface-size(bedge,1),1);
                v2=inedge(iface-size(bedge,1),2);
                IJ=coord(v2,:)-coord(v1,:);
                norma=norm(IJ);
                nij=R*IJ'/norma;
                auxpoint=(coord(v2,:)+coord(v1,:))*0.5;
                x=auxpoint(1,1);
                y=auxpoint(1,2);
                
                phi1=y-alfa*(x-0.5)-0.475;
                phi2=phi1-0.05;
                theta=atand(alfa);
                if phi1<0
                    
                    %Initialize "R":
                    
                    k = [1 0; 0 1];
                    
                    %Buld "kmap" again
                    KK = [k(1,1) k(1,2) 0;k(2,1) k(2,2) 0; 0 0 0];
                    
                elseif phi1>0 && phi2<0
                    
                    %Initialize "R":
                    
                    k = [0.01 0; 0 0.01];
                    
                    %Buld "kmap" again
                    KK= [k(1,1) k(1,2) 0; k(2,1) k(2,2) 0; 0 0 0];
                    a=-KK*[100*alfa;-100;0];
                    
                elseif phi2>0
                    
                    %Initialize "R":
                    
                    k = [1 0; 0 1];
                    
                    %Buld "kmap" again
                    KK = [k(1,1) k(1,2) 0; k(2,1) k(2,2) 0; 0 0 0];
                    a=-KK*[alfa;-1;0];
                end
                
                
            end
            F(iface,1)=dot(a,nij');
            
        end
        vel=F;
        K=kmap;
    case 'gaowu9'
        % adaptado de Gao e Wu, 2014
        auxbedge1=bedge(:,1);
        auxbedge2=bedge(:,2);
        
        bedge(:,1)=auxbedge2;
        bedge(:,2)=auxbedge1;
        
        u=0;
        alfa=0.2;
        for i=1:size(centelem,1)
            x=centelem(i,1);
            y=centelem(i,2);
            phi1=y-alfa*(x-0.5)-0.475;
            phi2=phi1-0.05;
            theta=atand(alfa);
            if phi1<0
                elem(i,5)=i;
                %Initialize "R":
                R = zeros(2);
                k = [1 0; 0 0.1];
                %Fill "R"
                R(1,1) = cosd(theta);
                R(1,2) = sind(theta);
                R(2,1) = -R(1,2);
                R(2,2) = R(1,1);
                %Fill "k" turning the tensor
                A=inv(R);
                k = A*k*R;
                %Buld "kmap" again
                kmap(i,1:5) = [i k(1,1) k(1,2) k(2,1) k(2,2)];
                
            elseif phi1>0 && phi2<0
                elem(i,5)=i;
                %Initialize "R":
                R = zeros(2);
                k = [100 0; 0 10];
                %Fill "R"
                R(1,1) = cosd(theta);
                R(1,2) = sind(theta);
                R(2,1) = -R(1,2);
                R(2,2) = R(1,1);
                %Fill "k" turning the tensor
                A=inv(R);
                k = A*k*R;
                %Buld "kmap" again
                kmap(i,1:5) = [i k(1,1) k(1,2) k(2,1) k(2,2)];
                fonte(i,1)=0;
            elseif phi2>0
                elem(i,5)=i;
                %Initialize "R":
                R = zeros(2);
                k = [1 0; 0 0.1];
                %Fill "R"
                R(1,1) = cosd(theta);
                R(1,2) = sind(theta);
                R(2,1) = -R(1,2);
                R(2,2) = R(1,1);
                %Fill "k" turning the tensor
                A=inv(R);
                k = A*k*R;
                %Buld "kmap" again
                kmap(i,1:5) = [i k(1,1) k(1,2) k(2,1) k(2,2)];
                fonte(i,1)= 0;
            end
            
            u(i,1)=2-x-alfa*y;
            KK=kmap(i,2:5);
            r=[KK(1,1) KK(1,2);
                KK(1,3) KK(1,4)];
            normKmap(i,1)=norm(r);
        end
        
        for iface=1:size(bedge,1)+size(inedge,1)
            R=[0 -1 0; 1 0 0; 0 0 0];
            if iface< size(bedge,1) || iface==size(bedge,1)
                v1=bedge(iface,1);
                v2=bedge(iface,2);
                IJ=coord(v2,:)-coord(v1,:);
                norma=norm(IJ);
                nij=R*IJ'/norma;
                auxpoint=(coord(v2,:)+coord(v1,:))*0.5;
                x=auxpoint(1,1);
                y=auxpoint(1,2);
                phi1=y-alfa*(x-0.5)-0.475;
                phi2=phi1-0.05;
                theta=atand(alfa);
                if phi1<0
                    
                    %Initialize "R":
                    R = zeros(2);
                    k = [1 0; 0 0.1];
                    %Fill "R"
                    R(1,1) = cosd(theta);
                    R(1,2) = sind(theta);
                    R(2,1) = -R(1,2);
                    R(2,2) = R(1,1);
                    %Fill "k" turning the tensor
                    A=inv(R);
                    k = A*k*R;
                    %Buld "kmap" again
                    KK = [k(1,1) k(1,2) 0;k(2,1) k(2,2) 0; 0 0 0];
                    
                elseif phi1>0 && phi2<0
                    
                    %Initialize "R":
                    R = zeros(2);
                    k = [100 0; 0 10];
                    %Fill "R"
                    R(1,1) = cosd(theta);
                    R(1,2) = sind(theta);
                    R(2,1) = -R(1,2);
                    R(2,2) = R(1,1);
                    %Fill "k" turning the tensor
                    A=inv(R);
                    k = A*k*R;
                    %Buld "kmap" again
                    KK= [k(1,1) k(1,2) 0; k(2,1) k(2,2) 0; 0 0 0];
                    
                elseif phi2>0
                    
                    %Initialize "R":
                    R = zeros(2);
                    k = [1 0; 0 0.1];
                    %Fill "R"
                    R(1,1) = cosd(theta);
                    R(1,2) = sind(theta);
                    R(2,1) = -R(1,2);
                    R(2,2) = R(1,1);
                    %Fill "k" turning the tensor
                    A=inv(R);
                    k = A*k*R;
                    %Buld "kmap" again
                    KK = [k(1,1) k(1,2) 0; k(2,1) k(2,2) 0; 0 0 0];
                    
                end
                
                a=-KK*[-1;-0.2;0];
            else
                v1=inedge(iface-size(bedge,1),1);
                v2=inedge(iface-size(bedge,1),2);
                IJ=coord(v2,:)-coord(v1,:);
                norma=norm(IJ);
                nij=R*IJ'/norma;
                auxpoint=(coord(v2,:)+coord(v1,:))*0.5;
                x=auxpoint(1,1);
                y=auxpoint(1,2);
                
                phi1=y-alfa*(x-0.5)-0.475+1e-8;
                phi2=phi1-0.05;
                theta=atand(alfa);
                if phi1<0
                    
                    %Initialize "R":
                    R = zeros(2);
                    k = [1 0; 0 0.1];
                    %Fill "R"
                    R(1,1) = cosd(theta);
                    R(1,2) = sind(theta);
                    R(2,1) = -R(1,2);
                    R(2,2) = R(1,1);
                    %Fill "k" turning the tensor
                    A=inv(R);
                    k = A*k*R;
                    %Buld "kmap" again
                    KK = [k(1,1) k(1,2) 0; k(2,1) k(2,2) 0; 0 0 0];
                    
                elseif phi1>0 && phi2<0
                    
                    %Initialize "R":
                    R = zeros(2);
                    k = [100 0; 0 10];
                    %Fill "R"
                    R(1,1) = cosd(theta);
                    R(1,2) = sind(theta);
                    R(2,1) = -R(1,2);
                    R(2,2) = R(1,1);
                    %Fill "k" turning the tensor
                    A=inv(R);
                    k = A*k*R;
                    %Buld "kmap" again
                    KK= [k(1,1) k(1,2) 0;k(2,1) k(2,2) 0;0 0 0];
                    
                elseif phi2>0
                    
                    %Initialize "R":
                    R = zeros(2);
                    k = [1 0; 0 0.1];
                    %Fill "R"
                    R(1,1) = cosd(theta);
                    R(1,2) = sind(theta);
                    R(2,1) = -R(1,2);
                    R(2,2) = R(1,1);
                    %Fill "k" turning the tensor
                    A=inv(R);
                    k = A*k*R;
                    %Buld "kmap" again
                    KK = [k(1,1) k(1,2) 0; k(2,1) k(2,2) 0;0 0 0];
                end
                
                a=-KK*[-1;-0.2;0];
                
            end
            F(iface,1)=dot(a,nij');
            
        end
        vel=F;
        K=kmap;
    case 'nikitin'
        for i=1:size(centelem,1)
            if centelem(i,2)<=(0.51-centelem(i,1))
                elem(i,5)=i;
                
                K(i,1:5)=[i 505 495 495 505];
            elseif (0.51-centelem(i,1))<centelem(i,2) && centelem(i,2)<(0.99-centelem(i,1))
                elem(i,5)=i;
                
                K(i,1:5)= [i 1000 0 0 10];
            elseif (0.99-centelem(i,1))<=centelem(i,2) && centelem(i,2)<(1.49-centelem(i,1))
                elem(i,5)=i;
                
                K(i,1:5)=[i 10 0 0 1000];
            elseif (1.49-centelem(i,1))<=centelem(i,2)
                elem(i,5)=i;
                
                K(i,1:5)=[i 505 495 495 505];
            end
        end
    case 'nikitin1'
        for i=1:size(centelem,1)
            if centelem(i,2)<=(0.251-centelem(i,1))
                elem(i,5)=i;
                
                K(i,1:5)=[i 505 495 495 505];
            elseif (0.251-centelem(i,1))<centelem(i,2) && centelem(i,2)<(0.51-centelem(i,1))
                elem(i,5)=i;
                
                K(i,1:5)= [i 10 0 0 1000];
            elseif (0.51-centelem(i,1))<=centelem(i,2) && centelem(i,2)<(0.751-centelem(i,1))
                elem(i,5)=i;
                
                K(i,1:5)=[i 1000 0 0 10];
                
                
            elseif (0.751-centelem(i,1))<=centelem(i,2) && centelem(i,2)<(0.99-centelem(i,1))
                elem(i,5)=i;
                
                K(i,1:5)=[i 10 0 0 1000];
            elseif (0.99-centelem(i,1))<=centelem(i,2) && centelem(i,2)<(1.249-centelem(i,1))
                elem(i,5)=i;
                
                K(i,1:5)=[i 1000 0 0 10];
            elseif (1.249-centelem(i,1))<=centelem(i,2) && centelem(i,2)<(1.49-centelem(i,1))
                elem(i,5)=i;
                
                K(i,1:5)=[i 10 0 0 1000];
            elseif (1.49-centelem(i,1))<=centelem(i,2) && centelem(i,2)<(1.749-centelem(i,1))
                elem(i,5)=i;
                
                K(i,1:5)=[i 1000 0 0 10];
            elseif (1.749-centelem(i,1))<=centelem(i,2)
                elem(i,5)=i;
                
                K(i,1:5)=[i 505 495 495 505];
                
            end
        end
    case 'lamine'
        K(1,1:5)=[1 5.5 4.5 4.5 5.5];
        elem(:,5)=1;
    case 'durlofsky'
        
        K(1,1:5)=[1 1 0 0 1];
        elem(:,5)=1;
    case 'shuec1'
        for i = 1:size(centelem,1)
            epsilon= rand(1,1);
            s= 0.1+ 50*(1+sin(10*(centelem(i,1)+centelem(i,2))))*epsilon;
            kmap(i,1:5)=[i s 0 0 s];
            K(i)=norm([kmap(i,2:3); kmap(i,4:5)]);
            
        end
    case 'shuec2'
        v=randperm(length(centelem));
        for i = 1:size(centelem,1)
            if centelem(i,1)<0.5
                Sum=0;
                for m=1:40
                    xl=centelem(v(m),:);
                    Sum=Sum + exp(-(norm(centelem(i,:)-xl)/0.05)^2);
                end
                
                kmap(i,1:5)=[i min(max(Sum,0.05),0.95) 0 0 min(max(Sum,0.05),0.95)];
                K(i)=norm([kmap(i,2:3); kmap(i,4:5)]);
                
            else
                
                epsilon= rand(1,1);
                s= 0.1+ 50*(1+sin(10*(centelem(i,1)+centelem(i,2))))*epsilon;
                kmap(i,1:5)=[i s 0 0 s];
                K(i)=norm([kmap(i,2:3); kmap(i,4:5)]);
                
            end
        end
    case 'shuec3'
        v=randperm(length(centelem));
        for i = 1:size(centelem,1)
            if centelem(i,1)<0.5 && centelem(i,2)<0.5
                epsilon= rand(1,1);
                s= 0.1+ 50*(1+sin(10*(centelem(i,1)+centelem(i,2))))*epsilon;
                kmap(i,1:5)=[i s 0 0 s];
                K(i)=norm([kmap(i,2:3); kmap(i,4:5)]);
                
            elseif centelem(i,1)>0.5 && centelem(i,2)<0.5
                
                Sum=0;
                for m=1:40
                    xl=centelem(v(m),:);
                    Sum=Sum + exp(-(norm(centelem(i,:)-xl)/0.05)^2);
                end
                
                kmap(i,1:5)=[i min(max(Sum,0.05),0.95) 0 0 min(max(Sum,0.05),0.95)];
                K(i)=norm([kmap(i,2:3); kmap(i,4:5)]);
                
            elseif centelem(i,1)<0.5 && centelem(i,2)>0.5
                
                Sum=0;
                for m=1:40
                    xl=centelem(v(m),:);
                    Sum=Sum + exp(-(norm(centelem(i,:)-xl)/0.05)^2);
                end
                
                kmap(i,1:5)=[i min(max(Sum,0.05),0.95) 0 0 min(max(Sum,0.05),0.95)];
                K(i)=norm([kmap(i,2:3); kmap(i,4:5)]);
                
            else
                
                epsilon= rand(1,1);
                s= 0.1+ 50*(1+sin(10*(centelem(i,1)+centelem(i,2))))*epsilon;
                kmap(i,1:5)=[i s 0 0 s];
                K(i)=norm([kmap(i,2:3); kmap(i,4:5)]);
                
            end
            
        end
    case 'altamenteheterogeneo'
        % mexer m e l
        for i=1:size(centelem,1)
            x=centelem(i,1);
            y=centelem(i,2);
            k0=1;
            s=4;
            m=3;
            l=3;
            r=k0*exp(sqrt(s)*cos(2*pi*m*x)*cos(2*pi*l*y));
            K(i,1:5) = [i r 0 0 r];
            normKmap(i,1)= r;
            elem(i,5)=i;
            u=0;
        end
    case 'pinchout'
        % não de definimos porque tudo esta definido no start
        elem=elem;
        K=kmap;
end
kmap=K;
fonte=fonte;
end

function teta=calculoteta(x,y)

if y<0
    if x>0
        m=atan(y/x);
    elseif x<0 && y>=0
        m=atan(y/x)+pi;
    elseif x<0 && y<0
        m=atan(y/x)-pi;
    elseif x==0 && y>0
        m=pi/2;
    elseif x==0 && y<0
        m=-pi/2;
    else
        m='indefined';
    end
    teta=2*pi+m;
else
    if x>0
        m=atan(y/x);
    elseif x<0 && y>=0
        m=atan(y/x)+pi;
    elseif x<0 && y<0
        m=atan(y/x)-pi;
    elseif x==0 && y>0
        m=pi/2;
    elseif x==0 && y<0
        m=-pi/2;
    else
        m='indefined';
    end
    teta=m;
    
end
end

function dtetadx=calculotetadx(x,y)

if x>0
    m=(-y)/(x^2+y^2);
elseif x<0 && y>=0
    m=(-y)/(x^2+y^2);
elseif x<0 && y<0
    m=(-y)/(x^2+y^2);
elseif x==0 && y>0
    m=0;
elseif x==0 && y<0
    m=0;
else
    m='indefined';
end
dtetadx=m;
end

function dtetady=calculotetady(x,y)

if x>0
    m=x/(x^2+y^2);
elseif x<0 && y>=0
    m=x/(x^2+y^2);
elseif x<0 && y<0
    m=x/(x^2+y^2);
elseif x==0 && y>0
    m=0;
elseif x==0 && y<0
    m=0;
else
    m='indefined';
end
dtetady=m;
end