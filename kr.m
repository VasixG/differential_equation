fimplicit3(@(x,y,z)(x+z)*y/(x*x)-1, [-6 6 -6 6 -6 6]);

hold on;

fimplicit3(@(x,y,z) log(abs(x+z)) - log(abs(z/x)) - log(2), [-6 -0.001 -6 6 -6 6]);
fimplicit3(@(x,y,z) log(abs(x+z)) - log(abs(z/x)) - log(2), [0.001 6 -6 6 -6 6]);

hold off;
