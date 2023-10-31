function [mask,xcontour,ycontour] = blWeights(bf)
    drise = 1e-3; % Size of smooth step

    [nx,ny] = size(bf.X);
    vortz = bf.Vortz;
    X = bf.X;
    Y = bf.Y;
    
    xcontour = zeros(nx,1);
    ycontour = zeros(nx,1);
    mask = zeros(size(X));
    for i=1:nx
        xi = X(i,:);
        yi = Y(i,:);
        di = sqrt((xi-xi(1)).^2+(yi-yi(1)).^2);
        vi = abs(vortz(i,:));
        vortThr = 0.03*vi(1); % Vorticity threashold
    
        for j=1:ny-1
            if vortThr < vi(j) && vortThr > vi(j+1)
                xcontour(i) = interp1([vi(j),vi(j+1)],[xi(j),xi(j+1)],vortThr);
                ycontour(i) = interp1([vi(j),vi(j+1)],[yi(j),yi(j+1)],vortThr);
                break
            end
        end
     
        dT = sqrt((xcontour(i)-xi(1))^2+(ycontour(i)-yi(1))^2);
        mask(i,:) = 1-step(di,dT-drise/2,drise);
    end
end

function y = step(xx,xstart,drise)
    y = zeros(size(xx));
    idx_after = find(xx>=(xstart+drise));
    idx = find(xx>xstart & xx<(xstart+drise));
    y(idx_after) = ones(length(idx_after),1);
    xnew = (xx(idx)-xstart)/drise;
    y(idx) = 1./(1+exp(1./(xnew-1)+1./xnew));
end
