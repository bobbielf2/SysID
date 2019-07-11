figure, set(gcf,'position',[200,200,600,300],'color','w')
for n = 2:7
    subplot(2,3,n-1)
    d = 0.4; s = 0.7;
    g = @(u,v) 3/8*((2*abs(u+1i*v).^(floor(n/2)*2)-1).^2-1) - imag(((u+1i*v)).^n);
    [X, Y] = meshgrid(linspace(-1,1,1000));
    Z = d*g(X/s,Y/s); Z(Z>0.05)=NaN;
    surf(X,Y,Z,'edgealpha',0)
    zmin = min(Z,[],'all');
    zmax = max(Z,[],'all');
    daspect([2.5,2.5,1])
    pbaspect([2,2,1])
    axis tight
    % zlim([zmin - (zmax-zmin)*0.2, 0])
    zlim([zmin, 0])
    title([num2str(n),'-well potential'])
    camlight
    camlight
end
%export_fig fig.pdf -opengl