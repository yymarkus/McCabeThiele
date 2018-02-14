function [xtrans, NTU] = PackedDistTransition(xf, xd, xb_1, R, F, D, B, run_number, Saving, NTU, q)

if NTU == true
%% Number of Transfer Units
xtrans_oldest = 0.11;
xtrans_old = 0.1;
%xb_1 = xb(xf, xd, F, D);
enr_oldest = NTUenr(xd, xtrans_oldest, R);
str_oldest = NTUstr(xd, xf, xtrans_oldest, xb_1, R, F, D);
enr_old = NTUenr(xd, xtrans_old, R);
str_old = NTUstr(xd, xf, xtrans_old, xb_1, R, F, D);
deltaNTU_old = str_old - enr_old;
deltaNTU_oldest = str_oldest - enr_oldest;
iterations = 0;
data_array = zeros(100, 3);
xtrans_new = xtrans_old;
while abs(deltaNTU_old) > 0.01 
    iterations = iterations + 1;
    xtrans_new = xtrans_old - (xtrans_old - xtrans_oldest)*deltaNTU_old/(deltaNTU_old - deltaNTU_oldest);
    xtrans_oldest = xtrans_old;
    xtrans_old = xtrans_new;
    enr_oldest = enr_old;
    str_oldest = str_old;
    enr_old = NTUenr(xd, xtrans_old, R);
    str_old = NTUstr(xd, xf, xtrans_old, xb_1, R, F, D);
    deltaNTU_oldest = str_oldest - enr_oldest;
    deltaNTU_old = str_old - enr_old;
    data_array(iterations, 1) = iterations;
    data_array(iterations, 2) = str_old;    
    data_array(iterations, 3) = enr_old;

end
figure(run_number*2-1);
plot (data_array(1:iterations,1), data_array(1:iterations,2), '-bs');
hold on
plot (data_array(1:iterations,1), data_array(1:iterations,3), '-rs');
titletext1 = ['Iterations of NTU Reflux Ratio ', num2str(R), ':1'];
titletext2 = ['Feed mol%: ' , num2str(round(xf*100,2)), ', Bottoms mol%: ', num2str(round(xb_1*100,2)), ', Distillate mol%: ', num2str(round(xd*100, 2))];
title({titletext1; titletext2});
legend('Stripping', 'Enriching')
xlabel('Iteration #')
ylabel('Number of Transfer Units')
if Saving == true
    saveas(figure(run_number*2-1), ['NTU Reflux: ', R, ' ', titletext2, '.png']);
    hold off
    close(figure(run_number*2-1));
end
end
%% McCabe Thiele
array_len = 1000;
x1 = zeros(array_len,1);
y1 = zeros(array_len,1);
y2 = zeros(array_len,1);
x1 = linspace(0, 1, array_len);
line_str = zeros(array_len,1);
line_enr = zeros(array_len,1);
for i=1:array_len
    y1(i) = fstar(x1(i));
    line_str(i) = yStripping(x1(i), xf, xd, R, F, D, xb_1, B);
    line_enr(i) = yEnriching(x1(i), xd, R);
end

figure(run_number*2);
x1 = x1';
p1 = plot(x1, y1); %Equilibrium plot
hold on
titletext1 = ['x-y Diagram Reflux Ratio ', num2str(R), ':1'];
titletext2 = ['Feed mol%: ' , num2str(round(xf*100,2)), ', Bottoms mol%: ', num2str(round(xb_1*100,2)), ', Distillate mol%: ', num2str(round(xd*100, 2))];
title({titletext1; titletext2});

McCabeIntersect_fun = @(x) [yStripping(x,xf,xd,R,F,D,xb_1, B) - yEnriching(x, xd, R)];
McCabeIntersect = fzero(McCabeIntersect_fun, xf);

if q ~= 1
    str_lower_limit = find(x1 >= xb_1, 1)-1;
    str_upper_limit = find(y1 > yStripping(McCabeIntersect,xf,xd,R,F,D,xb_1,B)+0.05, 1);

    p2 = plot(x1(str_lower_limit:str_upper_limit, :), line_str(str_lower_limit:str_upper_limit, :), 'Color',[0.5 0.2 0.8]);
end
axis([0 1 0 1]);
xlabel('x [Liquid mol% Ethanol]');
ylabel('y [Vapour mol% Ethanol]');


j=1;
xp = zeros(2,1);
yp = zeros(2,1);
xp(1)=xd;
yp(1)=xd;
y=xd;

if q ~= 1

enr_lower_limit = find(y1 > yEnriching(McCabeIntersect, xd, R), 1);
enr_upper_limit = find(line_enr >= xd, 1);
p3 = plot(x1(enr_lower_limit:enr_upper_limit), line_enr(enr_lower_limit:enr_upper_limit), 'Color', [0.9, 0.1, 0.6]);
set(line([0  1],[0  1]),'Color',[.5 .5 .5]);

while (xp(j)>McCabeIntersect),
%while (xp(j) > xf),
    fstar_fun = @(x) [fstar(x)-yp(j)];
    xp(j+1) = fzero(fstar_fun, xp(j));
    yp(j+1)= yEnriching(xp(j+1),xd,R);
    set(line([xp(j) xp(j+1)],[yp(j) yp(j)]),'Color',[0 0 1]);
    if (xp(j+1)>McCabeIntersect) set(line([xp(j+1) xp(j+1)],[yp(j) yp(j+1)]),'Color',[0 0 1]);
    end
        j=j+1;
end    
if j == 1
    yp(1) = xd;
    j = j+1;
end
yp(j) = yp(j-1);
while (xp(j)>xb_1),
    yp(j+1)= yStripping(xp(j),xf, xd, R, F, D, xb_1, B);
    fstar_fun = @(x) [fstar(x)-yp(j+1)];
    xp(j+1) = fzero(fstar_fun, 0.1);
    set(line([xp(j) xp(j)],[yp(j) yp(j+1)]),'Color',[0 0 1]);
    set(line([xp(j) xp(j+1)],[yp(j+1) yp(j+1)]),'Color',[0 0 1]);
    
    j=j+1;
end    
yp(j+1)= yStripping(xp(j),xf, xd, R, F, D, xb_1, B);
fstar_fun = @(x) [fstar(x)-yp(j+1)];
xp(j+1) = fzero(fstar_fun, 0.1);
set(line([xp(j) xp(j)],[yp(j) yp(j+1)]),'Color',[0 0 1]);
legend([p1 p2 p3],{'Equilibrium','Stripping','Enriching'});

qslope = (yStripping(McCabeIntersect, xf, xd, R, F, D, xb_1, B)-xf)/(McCabeIntersect-xf)


else
    
    McCabeIntersect = xf;
    enr_lower_limit = find(y1 > yEnriching(McCabeIntersect, xd, R), 1);
    enr_upper_limit = find(line_enr >= xd, 1);
    p3 = plot(x1(enr_lower_limit:enr_upper_limit), line_enr(enr_lower_limit:enr_upper_limit), 'Color', [0.9, 0.1, 0.6]);
    set(line([0  1],[0  1]),'Color',[.5 .5 .5]);
    
    while (xp(j)>McCabeIntersect),
    %while (xp(j) > xf),
        fstar_fun = @(x) [fstar(x)-yp(j)];
        xp(j+1) = fzero(fstar_fun, xp(j));
        yp(j+1)= yEnriching(xp(j+1),xd,R);
        set(line([xp(j) xp(j+1)],[yp(j) yp(j)]),'Color',[0 0 1]);
        if (xp(j+1)>McCabeIntersect) set(line([xp(j+1) xp(j+1)],[yp(j) yp(j+1)]),'Color',[0 0 1]);
        end
            j=j+1;
    end    
    if j == 1
        yp(1) = xd;
        j = j+1;
    end
    yp(j) = yp(j-1);
    slope = (yEnriching(xf, xd, R) - xb_1)/(xf-xb_1);
    y_intercept = xb_1-xb_1*slope;
    strip_fun = @(x) [slope.*x + y_intercept];
    x_4 = linspace(xb_1, xf);
    y_4 = strip_fun(x_4);
    p2 = plot(x_4, y_4, 'Color',[0.5 0.2 0.8]);
    legend([p1 p2 p3],{'Equilibrium','Stripping','Enriching'});
   
    
    while (xp(j)>xb_1),
        yp(j+1)= strip_fun(xp(j));%yStripping(xp(j),xf, xd, R, F, D, xb_1, B);
        fstar_fun = @(x) [fstar(x)-yp(j+1)];
        xp(j+1) = fzero(fstar_fun, 0.1);
        set(line([xp(j) xp(j)],[yp(j) yp(j+1)]),'Color',[0 0 1]);
        set(line([xp(j) xp(j+1)],[yp(j+1) yp(j+1)]),'Color',[0 0 1]);
    
        j=j+1;
    end    
    yp(j+1)= strip_fun(xp(j));%yStripping(xp(j),xf, xd, R, F, D, xb_1, B);
    fstar_fun = @(x) [fstar(x)-yp(j+1)];
    xp(j+1) = fzero(fstar_fun, 0.1);
    set(line([xp(j) xp(j)],[yp(j) yp(j+1)]),'Color',[0 0 1]);
    
    qslope = (strip_fun(McCabeIntersect)-xf)/(McCabeIntersect-xf)

end

q = (-1/(1/qslope -1))

if Saving == true
    saveas(figure(run_number*2), ['McCabe Reflux: ', R, ' ', titletext2, '.png']); 
    close(figure(run_number*2));
end
if NTU == true
    xtrans = xtrans_new;
    NTU = mean([str_old enr_old])*2;
else 
    xtrans = 0;
    NTU = 0;
end

end
