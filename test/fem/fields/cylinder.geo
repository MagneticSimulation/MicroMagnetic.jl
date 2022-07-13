#
## a cylinder
#
algebraic3d

# cut cylinder by planes:

solid fincyl = cylinder ( 0, 0, 5; 0, 0, 1; 2.0 )
	and plane (0, 0, 3; 0, 0, 1)
	and plane (0, 0, 1; 0, 0, -1) -maxh=0.4;

tlo fincyl;q
