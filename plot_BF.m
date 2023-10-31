clear all
data = base_case('fringe_m50.f00001',230,30);


XA = data.xx(1,:); % x coordinate points at the surface
YA = data.yy(1,:); % y coordinate points at the surface
[val0 inx0] = min(abs(XA));
dth = data.dth(ix0:end);
xa = XA(ix0:end);

figure()
plot(xa,dth)
xlabel('x')
ylabel('\delta')

