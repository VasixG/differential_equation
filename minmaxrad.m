c = 0.03;

f1 =@(x,y) x*y-1-0.03*exp((x^2 + y^2)/2);
f2 =@(x,y) x*y-1+0.03 * exp((x^2 + y^2)/2);

fimplicit(f1, [-5 5 -5 5]);

hold on;

fimplicit(f2, [-5 5 -5 5]);

fimplicit(@(x,y) x-y, [-5 5 -5 5], "DisplayName" ,"y = x");
fimplicit(@(x,y) x+y, [-5 5 -5 5], "DisplayName" ,"y = -x");

fimplicit(@(x,y) x^2+y^2 - (1.39)^2, [-4 4 -4 4], "--");

fimplicit(@(x,y) x^2+y^2 - (3.5)^2, [-4 4 -4 4], "--");

fimplicit(@(x,y) x^2+y^2 - (1.39)^2, [-4 4 -4 4], "--");

fimplicit(@(x,y) x^2+y^2 - (3.5)^2, [-4 4 -4 4], "--");

plot(sqrt(2), sqrt(2), ".");

plot(-sqrt(2), -sqrt(2), ".");

plot(0, 0, ".");
%, "x^2 + y^2 = (r_{*})^2", "x^2 + y^2 = (r^{*})^2"
lgd = legend("xy-1=3*10^{-2}e^{(x^2+y^2)/2}","xy-1=-3*10^{-2}e^{(x^2+y^2)/2}", "y = x", "y = -x");
lgd.Location = "southeastoutside";
hold off;
