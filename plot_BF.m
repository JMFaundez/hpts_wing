clear all
data = base_case('fringe_m50.f00001',230,30);


XA = data.xx(1,:); % x coordinate points at the surface
YA = data.yy(1,:); % y coordinate points at the surface
[val0 ix0] = min(abs(XA));
dth = data.dth(ix0:end);
d99 = data.y99(ix0:end);
xa = XA(ix0:end);
ya = YA(ix0:end);

figure()
hold on
plot(xa,dth)
plot(xa,d99)
xlabel('x')
ylabel('\delta')


X = data.xx;
Y = data.yy;

figure()
hold on
contourf(X,Y,data.ux,30,'LineStyle','none')
colorbar()
plot(xa,ya+dth)
axis equal
xlabel('x')
ylabel('y')


