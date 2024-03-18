%衍射函数
function DiffractZ=Dif(Object,lambda,dx,Zd)
%global lambda dx dy M N;
%Zc=N*dx^2/lambda;
delta=dx;
[M,N]=size(Object);
dy=dx;
Zc=M*dx^2/lambda;
%if abs(Zd)>Zc   
% %% 点扩散函数 
% [x,y] = meshgrid(-(ceil(M/2)-0.5)*delta:delta:(ceil(M/2)-0.5)*delta,-(ceil(N/2)-0.5)*delta:delta:(ceil(N/2)-0.5)*delta) ;%偶数使用这部分
% SphFunct=(1/(1i*lambda*Zd))*exp(1i*2.0*pi*Zd/lambda)*exp(1i*pi*(x.^2+y.^2)/(lambda*Zd));
% Object_F=fftshift(fft2(ifftshift(Object)));
% SphFunct_F=fftshift(fft2(ifftshift(SphFunct)));
% DiffractZ=(dx*dy).*fftshift(ifft2(ifftshift(SphFunct_F.*Object_F)));
%%
% else
%% 角谱法
du=1/(N*dx);dv=1/(M*dy);
[u,v]=meshgrid(-du*(floor(N/2)-0.5):du:du*(floor(N/2)-0.5),-dv*(floor(M/2)-0.5):dv:dv*(floor(M/2)-0.5)); %偶数
TransFunct=exp(1i*2.0*pi*Zd*((1/lambda)^2-u.^2-v.^2).^0.5);
Object_F=fftshift(fft2(ifftshift(Object)));
DiffractZ=fftshift(ifft2(ifftshift(Object_F.*TransFunct)));
%%
% end
end

