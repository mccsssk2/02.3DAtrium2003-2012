
clear all;
ss = sprintf('courtVW186_3.dat');
a = load(ss);
sizzee = size(a);
x=1:1:sizzee(2);
y=1:1:sizzee(1);
[X,Y]=meshgrid(x,y);
%sxy = sprintf('court2d%d',i);
surf(X,Y,a);
view(0,90);
caxis([-70 30]);
shading interp;
colormap(jet);
axis off;
grid off;
