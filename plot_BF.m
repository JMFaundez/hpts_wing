clear all
data = base_case('fringe_m50.f00001',230,30);


XA = data.xx(1,:); % x coordinate points at the surface
YA = data.yy(1,:); % y coordinate points at the surface
[val0 ix0] = min(abs(XA));
dth = data.dth(ix0:end);
d99 = data.y99(ix0:end);
xa = XA(ix0:end);
ya = YA(ix0:end);


bf_diego = load('baseFlow.mat');

[mask,xcontour,ycontour] = blWeights(bf_diego);


d_diego = sqrt((xcontour-bf_diego.X(:,1)).^2+(ycontour-bf_diego.Y(:,1)).^2) 

figure()
hold on
plot(xa,3*real(dth),'DisplayName','3\delta^*')
plot(xa,d99,'DisplayName','my d99')
plot(xcontour,d_diego,'DisplayName','diego')
legend()
xlabel('x')
ylabel('\delta')


X = data.xx;
Y = data.yy;

figure()
hold on
contourf(X,Y,data.ux,30,'LineStyle','none')
colorbar()
plot(xa,ya'+3*real(dth),'r--','LineWidth',1.5)
plot(xcontour,ycontour,'g-','LineWidth',1.5)
axis equal
xlabel('x')
ylabel('y')


