function p2 = script_roation( p1, finalaz, finalel)
    [az1,el1,r1] = cart2sph(p1(1),p1(2),p1(3));
    [x,y,z] =  sph2cart(finalaz,finalel,r1);
    [az1test,eltest,r1test] = cart2sph(x,y,z);
     p2 = [x,y,z];




