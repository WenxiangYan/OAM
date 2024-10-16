%% Global OAM for vortex-free self-accelerating beams with the unclosed parabolic-linear trajectory in Fig.3(U-X)
clear  
clc;
%% Basic parameter  
wave_length=532e-9;
obj_sampling=2048;
pixel_size=3.74e-6;
d=linspace(-(obj_sampling/2-0.5)*pixel_size,(obj_sampling/2-0.5)*pixel_size,obj_sampling); 
[X,Y]=meshgrid(d,-d);dx=X(1,2)-X(1,1);
[fai,r]=cart2pol(X,Y);
color_hot=colormap(hot);
mycolor=[color_hot(:,3),color_hot(:,1),color_hot(:,2)];
f=0.2; % focal length of Optical Fourier transform lens 
k0=2*pi/wave_length;kx=2*pi*X/(wave_length*f);ky=2*pi*Y/(wave_length*f);
[xita,kr]=cart2pol(kx,ky);
mask=CFO_circ_pixel(obj_sampling,obj_sampling,0,0,1024);
kz=sqrt(k0^2-kr.^2);kz(k0^2-kr.^2<=0)=0;
%% Theoretical predictions of local and global OAM by Eq.(5)
num=401; z1=linspace(-f,f,num); t=z1/f+1; 
S0=200e-6;ht=0.5*t;gt=-t.*(t-2); %the parabolic-linear trajectory in Fig.3(U-X)
gt1=gradient(gt,z1);ht1=gradient(ht,z1);
OAM_local_theory=k0*S0^2*(gt.*ht1-ht.*gt1);

%figure;plot(z1,OAM_local_theory);title('Theoretical prediction of local OAM');xlabel('z')
Iz=ones(1,num); % uniform intensity profile
OAM_global_theory=sum(OAM_local_theory.*ones(1,num))/num

%% Angular Specturm of Self-accelerating Bessel-like beams by Eq.(8)
kr0=0.25*max(max(kr))/sqrt(2); %the magnitude of the transverse component of the wavevector of the plane waves comprising the Bessel-like beam
scale=kr0/k0; kz0=k0*sqrt(1-scale^2);
z=linspace(-f,f,401); dz=z(2)-z(1);

E3=0;
for i=1:401 
t=z(i)/f+1; 
S0=200e-6;ht=0.5*t;gt=-t.*(t-2); % the parabolic-linear trajectory in Fig.3(U-X)
E=1;%uniform intensity profiles
E1=E.*exp(1i*(kz0)*z(i)).*exp(1i.*kx.*S0.*gt+1i.*ky.*S0.*ht); 
E2=E1.*exp(-1i*kz* z(i))*dz;
E3=E3+E2;
end
UKZ=E3./(2*pi*kz).*mask; % Angular Specturm
figure;subplot(1,2,1);imagesc((abs(UKZ).^2/max(max(abs(UKZ).^2))));title('Intensity of Angular Specturm');hold on
subplot(1,2,2);imagesc(angle(UKZ));title('Phase of Angular Specturm');
%% Direct calculation of global OAM from the angular spectrum by Eq.(9)
[Ukx,Uky]=gradient(UKZ,2*pi*d/(wave_length*f),-2*pi*d/(wave_length*f)); 
Ukf=kr.*(Uky.*cos(xita)-Ukx.*sin(xita)); 
Ck=-conj(UKZ).*1i.* Ukf; % OAM density
Energyk=abs(UKZ).^2; % intensity
OAM_global_calculation_from_angular_spectrum=real(sum((Ck(:)))/sum(Energyk(:))) 

%% Generation of Self-accelerating Bessel-like beams by Optical Fourier transform
b2=2048;UKZa=padarray(UKZ,[b2,b2]);
maskca=CFO_circ_pixel(obj_sampling+2*b2,obj_sampling+2*b2,0,0,2048);
UKZa=UKZa.*maskca;
obj_sampling1=obj_sampling+2*b2; %sampling for lens
d1=linspace(-dx*(obj_sampling1-1)/2,dx*(obj_sampling1-1)/2,obj_sampling1);
[X1,Y1]=meshgrid(d1,-d1); 
maskls=CFO_circ_pixel(obj_sampling1,obj_sampling1,0,0,2048); 
tj=exp(-1i*pi*(X1.^2+Y1.^2)/wave_length/f); % lens
DiffractZf=Dif(UKZa,wave_length,dx,f); % Dif:diffraction function based on angular spectrum theory
DiffractZb=DiffractZf.*tj.*maskls; 
zoom2=2;
Object_UKZ=DiffractZb(obj_sampling1/2-obj_sampling/zoom2+1:obj_sampling1/2+obj_sampling/zoom2,obj_sampling1/2-obj_sampling/zoom2+1:obj_sampling1/2+obj_sampling/zoom2);
% figure;imagesc(abs(Object_UKZ));colormap(mycolor);
% figure;imagesc(angle(Object_UKZ));colormap(jet);
nz_slice=101; z_distance_list=linspace(0,2*f,nz_slice);  % slice number along z between (-f,f)
DiffractZ3D=zeros(obj_sampling,obj_sampling,nz_slice); 
for ii=1:nz_slice 
Zd=z_distance_list(ii);DiffractZtemp=Dif(Object_UKZ,wave_length,dx,Zd);
DiffractZ3D(:,:,ii)=DiffractZtemp; 
end
Imax=max(max(max(abs(DiffractZ3D).^2)));
%% Observation of Self-accelerating Bessel-like beams
width=256;
for ii=1:nz_slice 
figure(200);
colormap(mycolor);
tt=abs(DiffractZ3D(1024-width+1:1024+width,1024-width+1:1024+width,ii).^2)/Imax; 
imagesc(tt);hold on
colorbar('FontSize',16);title(['Z=',num2str(z_distance_list(ii)*1000-200),'mm'],'FontSize',18);caxis([0 1]);
axis off
set(gcf,'unit','centimeters','position',[3 3 20 20]); 
set(gca,'Position',[.03 .03 0.83 0.83]); 
set(gcf,'color','w');

fig = getframe(gcf); 
imind = frame2im(fig);
[imind,cm] = rgb2ind(imind,256);
if ii == 1
name=' self-accelerating Bessel-like beam with parabolic-linear acceleration and uniform intensity profile.gif';
imwrite(imind,cm,name,'gif', 'LoopCount',inf,'DelayTime',0.1);
else
imwrite(imind,cm,name,'gif','WriteMode','append','DelayTime',0.1);
end
end
%% %% Computation of OAM using the complex spatial distribution in real space by Eq.(10) with allowances for errors from finite aperture dimensions
OAM_global_calculation_from_complex_spatial_distribution=zeros(1,nz_slice);
for ii=1:nz_slice
DiffractZd1=DiffractZ3D(:,:,ii)/sqrt(Imax);
[OAM_density,OAM_per_photon]=OAMz(DiffractZd1,zeros(obj_sampling,obj_sampling),zeros(obj_sampling,obj_sampling),d,-d);
OAM_global_calculation_from_complex_spatial_distribution(ii)=OAM_per_photon;
end
figure;plot(z_distance_list-f,OAM_global_calculation_from_complex_spatial_distribution,'r');
ylim([0 1.05*max(OAM_global_calculation_from_complex_spatial_distribution)]);
title('global OAM calculation from complex spatial distribution');xlabel('z')