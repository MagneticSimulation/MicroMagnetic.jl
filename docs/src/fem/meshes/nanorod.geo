algebraic3d

solid p1 = plane (0, 0, 1; 0, 0, -1) -bc=1;
solid p2 = plane (0, 0, 100; 0, 0, 1) -bc=2;

solid p3 = plane (0, 0, -1; 0, 0, 1) -bc=3;
solid p4 = plane (0, 0, -100; 0, 0, -1) -bc=4;

solid c1 = cylinder ( 0, 0, -100; 0, 0, 100; 5.0 )
	and p1 and p2 -maxh = 5;

solid c2 = cylinder ( 0, 0, -100; 0, 0, 100; 5.0 )
	and p3 and p4 -maxh = 5;

solid main = c1 or c2;
tlo main;

#identify periodic p1 p3;