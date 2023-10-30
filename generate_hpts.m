%close all
clear all

addpath matlab_script/

A = 0.2969;
B = 0.1260;
C = 0.3516;
D = 0.2843;
E = 0.1015;

naca = @(x) 5*0.08*(A*sqrt(x) - B*x - C*x.^2 + D*x.^3 - E*x.^4);
dnaca = @(x) 5*0.08*(0.5*A*1/sqrt(x) - B - 2*C*x + 3*D*x.^2 - 4*E*x.^3);

xf = 0.55;
xi = 0.01;
step = 8;

data = base_case('fringe_m50.f00001',230,30);
XA = data.xx(1,:); % x coordinate points at the surface
YA = data.yy(1,:); % y coordinate points at the surface
[val0 inx0] = min(abs(XA));
inxf = find(XA>=xf,1,'first');
XA = XA(inx0:inxf);
YA = YA(inx0:inxf);
inxi = find(XA>=xi,1,'first');
XA = XA(inxi:step:end);
YA = YA(inxi:step:end);

figure()
plot(XA,YA,'.')
axis('equal')

%x = linspace(0.002, 0.34,100);
%y = naca(x);
x = XA;
y = YA; 
[d,D1,D2,W] = chebmat_trans(70,15e-3,5e-3);
%d = linspace(0,0.08,5)*1e-3;
z = linspace(-0.035, 0.035, 61);
z = z(1:end-1);


nx = length(x);
nd = length(d);
nz = length(z);
Ntot = nx*nd*nz;

points = zeros(Ntot,3);

for i=1:nx
  iix = (i-1)*nd*nz+1;
  ifx = i*nd*nz;
  alpha = atan(dnaca(x(i)));
  for j=1:nd
	iiy = nz*(j-1);
	ify = nz*j-1;
	points(iix+iiy:iix+ify,1) = x(i)-d(j)*sin(alpha);
	points(iix+iiy:iix+ify,2) = y(i) + d(j)*cos(alpha);

	for k=1:nz
	  points(iix+iiy+k-1,3) = z(k);
	end
  end
end
fileID = fopen('FST.his','w');
fprintf(fileID,'%i\n', Ntot);
fprintf(fileID,'%1.15f %1.15f %2.8f\n',points');

%%
figure()
hold on
mesh(data.xx,data.yy,data.xx*0)
plot(points(1:nz:end,1),points(1:nz:end,2),'b.')
view(2)
axis equal

%%
figure()
plot3(points(:,3), points(:,1),points(:,2),'.')
axis('equal')




