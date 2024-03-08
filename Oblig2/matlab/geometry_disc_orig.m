% box of draught D and length (width) L, with D as unit length
% left vertical side
D=1;
L=2;
Nside=10;
Nbott=20;
N=Nside+Nbott+Nside;
dy=D/Nside;
dx=L/Nbott;
xp(1)=-L/2;
xm(1)=-L/2;
yp(1)=-dy;
ym(1)=-dy*(1-1);
coord=[xm(1),ym(1),xp(1),yp(1)];

for i=2:Nside;
xp(i)=-L/2;
xm(i)=-L/2;
yp(i)=-dy*i;
ym(i)=-dy*(i-1);
coord2=[xm(i),ym(i),xp(i),yp(i)];
coord=[coord;coord2];
end
for i=1+Nside:Nside+Nbott;
i1=i-Nside;
xp(i)=-L/2+dx*i1;
xm(i)=-L/2+dx*(i1-1);
yp(i)=-D;
ym(i)=-D;
coord2=[xm(i),ym(i),xp(i),yp(i)];
coord=[coord;coord2];
end
for i=1+Nside+Nbott:Nside+Nbott+Nside;
i1=i-Nside-Nbott;
xp(i)=L/2;
xm(i)=L/2;
yp(i)=-D+dy*i1;
ym(i)=-D+dy*(i1-1);
coord2=[xm(i),ym(i),xp(i),yp(i)];
coord=[coord;coord2];
end
save -ascii box.dat coord;
% plot geometry
xa=xm;
xa=[xm,xp(N)]
ya=ym;
ya=[ym,yp(N)]
hold on
axis ([-1.01 1.01 -1.1 0.01])
h1=plot(xa,ya,’k .’)
plot(xa,ya,’k -’)
axis off
get(h1);
set(h1,’MarkerSize’,[25])
set(gca,’FontSize’,20)