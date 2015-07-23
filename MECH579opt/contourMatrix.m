CONTZ=0;
minx=2.5;
maxx=4;
miny=3.2;
maxy=4;
divider=25;
minx=(divider*minx);
maxx=(divider*maxx);
miny=(divider*miny);
maxy=(divider*maxy);

for CONTI=-minx:maxx
    CONTI
    CONTX(CONTI+minx+1)=CONTI/(divider+1);
    for CONTJ=-miny:maxy
        CONTJ;
        points=[CONTI/(divider+1),CONTJ/(divider+1)];
        CONTZ(CONTI+minx+1,CONTJ+miny+1)=Contourf(points);
    end
end

CONTZ=CONTZ';
for CONTJ=-miny:maxy
    CONTY(CONTJ+miny+1)=CONTJ/(divider+1);
end

contour(CONTX,CONTY,CONTZ,75)