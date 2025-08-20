algebraic3d

solid sp = sphere(0, 0, 0; 10);

solid air = orthobrick(-2,-20,-20; 2, 20, 20) and not sphere(0, 0, 0; 12);

tlo sp -maxh=2.0;

tlo air -maxh=2.0;