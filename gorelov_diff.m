function gorelov_diff
format long;

x_0 = input("Enter x_0: ");
y_0 = input("Enter y_0: ");
t_0 = input("Enter t_0: ");

start_time = input("Enter start time: ");
end_time = input("Enter end time: ");

time_s = ceil(end_time-start_time)*1000;
root_s = 10000;
int_s = 10000;

time = linspace(start_time, end_time, time_s);

x = zeros(size(time));
y = zeros(size(time));

if x_0*y_0 == 1
    x_t = @(t) x_0*exp(t-t_0);
    y_t = @(t) 1./x_0*exp(-t+t_0);

    for i = 1:length(time)
        x(i) = x_t(time(i));
        y(i) = y_t(time(i));
    end
else
    
    if x_0 == 0 && y_0 == 0
        x_t = @(t) 0;
        y_t = @(t) 0;

        for i = 1:length(time)
            x(i) = x_t(time(i));
            y(i) = y_t(time(i));
        end
    end
    
    [theta_0,r_0] = cart2pol(x_0,y_0);

    c = 2*(x_0*y_0 - 1)*exp(-(x_0.^2+y_0.^2)/2);

    [R_ , Rp, Theta_, Thetap] = diffint(r_0, c, root_s, int_s, t_0, start_time, end_time, time_s);

    
    if theta_0 > -pi/4 && theta_0<= pi/4
        for i = 1:length(time)
            x(i) = Rp(i)*cos(Thetap(i));
            y(i) = Rp(i)*sin(Thetap(i));
        end
    end

    if theta_0 > pi/4 && theta_0<= 3*pi/4
        for i = 1:length(time)
            x(i) = R_(i)*cos(Theta_(i));
            y(i) = R_(i)*sin(Theta_(i));
        end
    end

    if (theta_0 > -pi && theta_0<= -3*pi/4) || (theta_0 > 3*pi/4 && theta_0< pi)
        for i = 1:length(time)
            x(i) = -Rp(i)*cos(Thetap(i));
            y(i) = -Rp(i)*sin(Thetap(i));
        end
    end

    if theta_0 > -3*pi/4 && theta_0 <= -pi/4
        for i = 1:length(time)
            x(i) = -R_(i)*cos(Theta_(i));
            y(i) = -R_(i)*sin(Theta_(i));
        end
    end

end

if x_0 == sqrt(2) && y_0 == sqrt(2)
    x_t = @(t) sqrt(2);
    y_t = @(t) sqrt(2);

    for i = 1:length(time)
        x(i) = x_t(time(i));
        y(i) = y_t(time(i));
    end
end

if x_0 == -sqrt(2) && y_0 == -sqrt(2)
    x_t = @(t) -sqrt(2);
    y_t = @(t) -sqrt(2);

    for i = 1:length(time)
        x(i) = x_t(time(i));
        y(i) = y_t(time(i));
    end
end

hold on;
xline(0);
yline(0);
fimplicit(@(x,y) x.*y-1 - c./2*exp((x.^2+y.^2)/2), [-10 10 -10 10]);
comet(real(x),real(y));
plot(x(1),y(1),"*", "MarkerSize",10);
% plot(time, R_);
% plot(time, Theta_);
hold off;
end

function [R_ , Rp, Teta_, Tetap] = diffint(r0, c, root_t, int_t, t0, time_start, time_end, time_t)

    f =@(x)1./sqrt(x.^2 - x.^(-2).*(c.*exp(x.^2/2)+2).^2);
    roots1 = @(r) (r^2-2)-c*exp(r^2/2);
    roots2 = @(r) (-r^2-2)-c*exp(r^2/2);

    h = @(r)r.^(-2) * (c*exp(r.^2/2) + 2);

    psi_p = @(r) -1/2 * asin(h(r)) + pi/2;
    psi_m = @(r) 1/2 *asin(h(r));

    list_ofroots1 = [];
    list_ofroots2re = [];
    list_ofroots1re = [];
    list_ofroots2 = [];

    t = linspace(time_start,time_end,time_t);

    R_ = zeros(size(t));
    Rp = zeros(size(t));

    Teta_ = zeros(size(t));
    Tetap = zeros(size(t));

    x = linspace(0,20, root_t);
    if c > 0
        for i = x
            val = fzero(roots1, i);
            if ~ ismember( round(val,4) , list_ofroots1) && ~isnan(val) && val >0
                list_ofroots1 = [list_ofroots1, round(val,4)];   
                list_ofroots1re = [list_ofroots1re, val];
            end
    
        end
    else
        for i = x
            val = fzero(roots1, i);
            if ~ ismember( round(val,4) , list_ofroots1) && ~isnan(val) && val > 0
                list_ofroots1 = [list_ofroots1, round(val,4)];   
                list_ofroots1re = [list_ofroots1re, val];
            end
    
            val2 = fzero(roots2, i);
            if ~ ismember( round(val2,4) , list_ofroots2) && ~isnan(val2) && val2 >0
                list_ofroots2 = [list_ofroots2, round(val2,4)];   
                list_ofroots2re = [list_ofroots2re, val2];
            end
    
        end
    end
    
    if c < 0
        if length(list_ofroots1re) == 1 && length(list_ofroots2re) == 1

            a = list_ofroots1re;
            b = list_ofroots2re;

            y = linspace(a,b, int_t);

            F = zeros(size(y));

            omega = real(integral(f,a,b));

            t0_ = real(integral(f,a,r0));
            t0p = real(integral(f,r0,b));

            for i = 1:length(y)
                F(i) = real(integral(f,r0,y(i)));
            end

            F_ = -F;

            for j = 1:length(t)
                k2 = floor((t(j)-t0_-t0)/omega);
                if mod(k2,2) == 0
                    R_(j) = interp1(F, y, t(j)-2*t0_-k2*omega-t0, "linear");
                    if mod(k2, 4) == 0
                        Teta_(j) = psi_m(R_(j));
                    end

                    if mod(k2, 4) == 2
                        Teta_(j) = psi_m(R_(j))+pi;
                    end
                else
                    R_(j) = interp1(F_, y, t(j)-(k2+1)*omega-t0, "linear");
                    if mod(k2, 4) == 1
                        Teta_(j) = psi_p(R_(j))+pi;
                    end
                    if mod(k2, 4) == 3
                        Teta_(j) = psi_p(R_(j));
                    end
                end
    
                k2p = floor((t(j)+t0_-t0)/omega);
    
                if mod(k2p,2) == 1
                    Rp(j) = interp1(F_, y, t(j)-2*t0p-(k2p-1)*omega - t0, "linear");
                    if mod(k2p, 4) == 1
                        Tetap(j) = psi_p(Rp(j))+pi;
                    end

                    if mod(k2p, 4) == 3
                        Tetap(j) = psi_p(Rp(j));
                    end
                else
                    Rp(j) = interp1(F, y, t(j)-(k2p)*omega -t0, "linear");
                    if mod(k2p, 4) == 0
                        Tetap(j) = psi_m(Rp(j));
                    end

                    if mod(k2p, 4) == 2
                        Tetap(j) = psi_m(Rp(j))+pi;
                    end
                end
    
            end
    
        end
    else
        if length(list_ofroots1re) == 2
            a = list_ofroots1re(1);
            b = list_ofroots1re(2);

            y = linspace(a,b, int_t);
            F = zeros(size(y));
            omega = real(integral(f,a,b));
            t0_ = real(integral(f,a,r0));
            t0p = real(integral(f,r0,b));
    
            for i = 1:length(y)
                F(i) = real(integral(f,r0,y(i)));
            end
    
            F_ = -F;
          
            for j = 1:length(t)
                k2 = floor((t(j)-t0_-t0)/omega);
                if mod(k2,2) == 0
                    R_(j) = interp1(F, y, t(j)-2*t0_-k2*omega-t0, "linear");
                    Teta_(j) = psi_m(R_(j));
                else
                    R_(j) = interp1(F_, y, t(j)-(k2+1)*omega-t0, "linear");
                    Teta_(j) = psi_p(R_(j));
                end
    
                k2p = floor((t(j)+t0_-t0)/omega);
    
                if mod(k2p,2) == 1
                    Rp(j) = interp1(F_, y, t(j)-2*t0p-(k2p-1)*omega-t0, "linear");
                    Tetap(j) = psi_p(Rp(j));
                else
                    Rp(j) = interp1(F, y, t(j)-(k2p)*omega-t0, "linear");
                    Tetap(j) = psi_m(Rp(j));
                end
    
            end
        end
    end
end
