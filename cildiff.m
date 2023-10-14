function cildiff
format long;

time_s = 10*1000;
root_s = 10000;
int_s = 10000;

time1 = linspace(0, 4.5, time_s);
time2 = linspace(0, 2.651439153110586, time_s);
time3 = linspace(0, 5.206071494613065, time_s);

x1 = zeros(size(time1));
y1 = zeros(size(time1));

x2 = zeros(size(time2));
y2 = zeros(size(time2));

x3 = zeros(size(time3));
y3 = zeros(size(time3));

c1 = 0;
c2 = 0.03;
c3 = -0.03;


[R_2 , Rp2, Teta_2, Tetap2] = diffint(0, c2*2, root_s, int_s, 0, 0, 2.651439153110586, time_s);
[R_3 , Rp3, Teta_3, Tetap3]= diffint(0, c3*2, root_s, int_s, 0, 0, 5.206071494613065, time_s);


x_t = @(t) 1/10*exp(t);
y_t = @(t) 10*exp(-t);

for i = 1:length(time1)
    x1(i) = x_t(time1(i));
    y1(i) = y_t(time1(i));
end

for i = 1:length(time2)
    x2(i) = -Rp2(i)*cos(Tetap2(i));
    y2(i) = -Rp2(i)*sin(Tetap2(i));
end

for i = 1:length(time3)
    x3(i) = -Rp3(i)*cos(Tetap3(i));
    y3(i) = -Rp3(i)*sin(Tetap3(i));
end
f1 = @(x,y,z) x.*(x>0).*y.*(y>0)-1 - c1./2*exp((x.^2+y.^2)/2);
f2 = @(x,y,z) x.*(x<0).*y.*(y<0)-1 - c2.*exp((x.^2.*(y<0)+y.^2.*(y<0))/2);
f3 = @(x,y,z) x.*y-1 - c3.*exp((x.^2+y.^2)/2);

%fimplicit3(f1, [-10 10 -10 10 -10 10],'FaceAlpha',.5, 'FaceColor','#f7d679', "EdgeColor","#94541c","LineWidth",0.1);
%hold on;

fimplicit3(f3, [-10 10 -10 10 0 7],'FaceAlpha',.2, 'FaceColor','#22311d', "EdgeColor","#22311d","LineWidth",0.2);
hold on;
fimplicit3(f2, [-10 10 -10 10 0 7],'FaceAlpha',.2, 'FaceColor','#0000b3', "EdgeColor","#0000b3","LineWidth",0.2);
fimplicit3(f1, [-10 10 -10 10 -1 10],'FaceAlpha',.05, 'FaceColor','#FF00FF', "EdgeColor","#94541c","LineWidth",0.3);
plot3(x3,y3,time3,"LineWidth",3, 'Color',"#22311d");
plot3(x2,y2,time2,"LineWidth",3, 'Color',"#0000b3");
plot3(x1,y1,time1,"LineWidth",3, 'Color',"#FF00FF");
% xline(0);
% yline(0);
% fimplicit(@(x,y) x.*y-1, [-2 2 -2 2]);
% plot(x1,y1);
% plot(x1(1),y1(1),"*", "MarkerSize",10);
% plot(time, R_);
% plot(time, Theta_);
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
            r0 = a;
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

            r0 = a;

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
