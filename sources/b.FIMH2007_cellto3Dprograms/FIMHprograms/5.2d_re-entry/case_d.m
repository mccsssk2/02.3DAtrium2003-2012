start=1;
finish=1000;
increment=5;
j = 1;
for i=start:increment:finish
    ss=sprintf('AF2/court_cai2d%d.dat',i);
    a=load(ss); 
    sizzee = size(a);
    x=1:1:sizzee(2);
    y=1:1:sizzee(1);
    [X,Y]=meshgrid(x,y);
    surf(X,Y,a);
    view(0,90);
    shading interp;
   % colormap(jet);
    caxis([0 0.0005]);
    axis([0 sizzee(2) 0 sizzee(1)]);
    M(:,j)=getframe;
    hold off;
    j = j+1;
end;
movie(M);
movie2avi(M,'af1_0.33.avi','fps',5,'quality',100);