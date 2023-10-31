% This script is used to create the interpolation mesh for the data.
clear all
%% Loading the data:

%coords = load('/scratch/kenzo/2021_FST_SPOD/for_BRAZIL/coordinates.mat');
%X = coords.output_x(:,:,1);
%Y = coords.output_y(:,:,1);
%Z = squeeze(coords.output_z(1,1,:));
coords = load("surface_points.mat");
%0.0418:8.4000e-04*5:0.3347;
dt_sn = 1.2e-3;
dx_p = dt_sn*5*0.7;

fact_x = 6/(1.5111e-5/1.204);
xf = 1.6e5/fact_x;
x0 = 0.2e5/fact_x;
step = 4;
%Xs = X(:,end,1);
%Ys = Y(:,end,1);
Xs = coords.xun;
Ys = coords.yun;
figure(1)
plot(Xs,Ys,'o')
axis equal

[val0, inx0] = min(abs(Xs));
ix = find(Xs(inx0:end)>=x0,1,'first');
inx0 = inx0+ix-10;
inxf = find(Xs>=xf,1,'first');

%x = linspace(x0,xf,40);
x = x0:dx_p:xf;
y = interp1(Xs(inx0:end),Ys(inx0:end),x,'nearest');
%x = Xs(inx0:step:inxf);
%y = Ys(inx0:step:inxf);

figure(1)
hold on
plot(x,y,'x')
axis equal
N=44;
ymax=0.01;
%[d,D1,D2,W] = chebmat_trans(100,0.15,0.1);
%[d,D1,D2,W] = chebmat(100,0.15);
beta = linspace(0,pi/2,N);
d = (1-cos(beta))/2*ymax*2;

%ybl = 0.002;
%d1 = linspace(0,ybl,N/2+1);
%d1 = (1-cos(beta))/2*ybl*2;
%d2 = linspace(0,1,N/2).^2*ymax;
%d2 = (1-cos(beta))/2*ymax*2;
%d2 = d2(2:end)+ybl;
%d = [d1,d2];

%d = linspace(0,1,N).^4*ymax;
figure()
hold on
plot(d,'k.')
%plot(d1,d1,'r.')
%plot(d2,d2,'bo')

z = linspace(-0.075, 0.075, 501);
z = z(1:end-1);
%z = 0;
nx = length(x);
nd = length(d);
nz = length(z);

Ntot = nx*nd*nz;
%%
%x1 =  X(:,end,1);
%x2 =  X(:,end-1,1);
%y1 =  Y(:,end,1);
%y2 =  Y(:,end-1,1);

% dist = sqrt( (x2-x1).^2 + (y2-y1).^2 );
% sin_a = (x2-x1)./dist;
% cos_a = (y2-y1)./dist;
% cos_a = cos_a(inx0:step:inxf);
% sin_a = sin_a(inx0:step:inxf);
% cos_a(1) = 0;
% sin_a(1) = -1;

%sin_a = coords.normal(inx0:step:inxf,1);
%cos_a = coords.normal(inx0:step:inxf,2);

cos_a = interp1(Xs(inx0:end),coords.normal(inx0:end,2),x);
sin_a = interp1(Xs(inx0:end),coords.normal(inx0:end,1),x);
cos_a=zeros(size(cos_a))+1;
sin_a=zeros(size(sin_a));

% figure()
% hold on
% plot(x,cos_a)
% plot(x,sin_a)
%plot(x,smooth(cos_a))

%%
points = zeros(Ntot,3);

for i=1:nx
  iix = (i-1)*nd*nz+1;
  ifx = i*nd*nz;
  for j=1:nd
	iiy = nz*(j-1);
	ify = nz*j-1;
	points(iix+iiy:iix+ify,1) = x(i)+d(j)*sin_a(i);
	points(iix+iiy:iix+ify,2) = y(i) + d(j)*cos_a(i);

	for k=1:nz
	  points(iix+iiy+k-1,3) = z(k);
	end
  end
end
fileID = fopen('sec_inst_mesh_V4.his','w');
fprintf(fileID,'%i\n', Ntot);
fprintf(fileID,'%1.15f %1.15f %2.8f\n',points');

%%
figure()
hold on
plot(points(1:nz:end,1),points(1:nz:end,2),'b.')
view(2)
axis equal

fact_x = 6/(1.5111e-5/1.204);
Rex0 = x0*fact_x
Rexf = xf*fact_x





