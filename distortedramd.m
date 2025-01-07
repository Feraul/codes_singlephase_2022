function [auxcoord]=distortedramd
global coord bedge

auxcoord=zeros(size(coord,1),4);

a=-1;
b=1;
% para malha quadilateral da SPE usei limites 10 e -10
% com alpha 0.8 e h=1/220
ex = a + (b-a)*rand(size(coord,1),1);
ey = a + (b-a)*rand(size(coord,1),1);
% quando mas proximo a 1 mas ditorcido a malha
alpha=0.3;
h=3; % coloque o valor do delta x= comprimento/numero de espacamentos 
%h=0.25 --> 160x80
%h=0.5 --> 80x40
%h=0.95 -->40x20
%h=2 -->20x10

for icoord=1:size(coord,1)
    x=coord(icoord,1);

    y=coord(icoord,2);
    if icoord>size(bedge,1) %&& abs(coord(icoord,1)-0.5)>1e-10
        % if single(y)~=0.5
        %   if icoord~=433 & icoord~=434 & icoord~=459 & icoord~=460
        %if icoord~=364 & icoord~=365 & icoord~=496 & icoord~=397 &...
        %        icoord~=380 & icoord~=381 & icoord~=412 & icoord~=413 &...
        %       icoord~=876 & icoord~=877 & icoord~=908 & icoord~=909 &...
        %        icoord~=892 & icoord~=893 & icoord~=924 & icoord~=925
           if (single(x)~=150) && (single(x)~=270) 
            x=coord(icoord,1)+h*alpha*ex(icoord,1);

            y=coord(icoord,2)+h*alpha*ey(icoord,1);

            auxcoord(icoord,1)=icoord;
            auxcoord(icoord,2)=x;
            auxcoord(icoord,3)=y;
        else
            auxcoord(icoord,1)=icoord;
            auxcoord(icoord,2)=x;
            auxcoord(icoord,3)=y;
        end
    else

        auxcoord(icoord,1)=icoord;
        auxcoord(icoord,2)=x;
        auxcoord(icoord,3)=y;

    end
end

end
