function [ Loz, SUM_Loz ] = OAMz( Ex, Ey, Ez, x_image, y_image )
%OAMz Orbital angular momentum, z component
if nargin<4
    x_image=1;
    y_image=1;
end
[X,Y]=meshgrid(x_image,y_image);
[Finxx,Finxy]=gradient(Ex,x_image,y_image);
[Finyx,Finyy]=gradient(Ey,x_image,y_image);
[Finzx,Finzy]=gradient(Ez,x_image,y_image);
Loz=imag(conj(Ex).*(X.*Finxy-Y.*Finxx)+conj(Ey).*(X.*Finyy-Y.*Finyx)+conj(Ez).*(X.*Finzy-Y.*Finzx));

Energy=sum(abs(Ex(:)).^2)+sum(abs(Ey(:)).^2)+sum(abs(Ez(:)).^2);
SUM_Loz=sum(sum(Loz))/Energy;
end

