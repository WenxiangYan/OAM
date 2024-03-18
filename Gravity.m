function [xm,ym] = Gravity(I3,n)
%输入的值要平方
Itotalx=0;
Itotaly=0;
I3=I3.^n;
Itotal=sum(sum(I3));
s=size(I3);
a=s(1);b=s(2);

[Y3,X3]=meshgrid(1:a,1:b);
Itotalx=sum(sum(X3.*I3));
Itotaly=sum(sum(Y3.*I3));

xm=Itotalx/Itotal;
ym=Itotaly/Itotal;
end