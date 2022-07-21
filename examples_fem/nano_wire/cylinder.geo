#
## a cylinder
#
algebraic3d

# cut cylinder by planes:

solid fincyl = cylinder ( 30, 0, 0; -1, 0, 0; 3 )
	and plane (0, 0, 0; -1, 0, 0)
	and plane (20, 0, 0; 1, 0, 0);

tlo fincyl -maxh=1.0;