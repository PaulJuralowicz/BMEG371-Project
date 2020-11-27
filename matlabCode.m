%% Paul's lovely numerical solution. We are working in micro meters and mili seconds
% The goal of this is to run fast, all actually analysis will be in the
% next block.

L = 20/1000;
h = 1 * 10^(-3);
k = 1 * 10^(-8);
N = cast((L/h)+1, 'int32');
x = 0:h:L;
D = 0.2;
Tmax = 0.1;
Tcurr = 0;
iterations = cast(Tmax/k,'int32'); %we want to goto 0.5 ms
w_prev = zeros(N,1);
w_curr = zeros(N,1);
Dkh = (D*k)/(h^2);
%initial condition
w_prev(1,1) = 1;
%BEGIN
t0 = clock;
for j = 1:iterations
    w_curr(1) = w_prev(1) + Dkh*(2*w_prev(2) - 2*w_prev(1));
    for i = 2:(N-1)
        w_curr(i) = w_prev(i) + Dkh*(w_prev(i-1) - 2*w_prev(i) + w_prev(i+1));
    end
    w_prev = w_curr;
end
ms = round(etime(clock,t0) * 1000);
%% Stuff to get nice graph
L = 20/1000;
h = 1 * 10^(-3);
k = 1 * 10^(-8);
N = cast((L/h)+1, 'int32');
x = 0:h:L;
D = 0.2;
Tmax = 0.1;
iterations = cast(Tmax/k,'int32'); %we want to goto 0.5 ms
w_prev = zeros(N,1);
w_curr = zeros(N,1);
Dkh = (D*k)/(h^2);
%initial condition
w_prev(1,1) = 1;
%BEGIN
weGotEm = false;
for j = 1:iterations
    w_curr(1) = w_prev(1) + Dkh*(2*w_prev(2) - 2*w_prev(1));
    for i = 2:(N-1)
        w_curr(i) = w_prev(i) + Dkh*(w_prev(i-1) - 2*w_prev(i) + w_prev(i+1));
    end
    w_prev = w_curr;
    if(abs(cast(j, 'like', k) * k - 0.00001) <0.00000002)
        w0001 = w_prev;
    end
    if(abs(cast(j, 'like', k) * k - 0.0003) <0.00000002)
        w0003 = w_prev;
    end
    if(abs(cast(j, 'like', k) * k - 0.0006) <0.00000002)
        w0006 = w_prev;
    end
    if(abs(cast(j, 'like', k) * k - 0.0009) <0.00000002)
        w0009 = w_prev;
    end
    if(abs(cast(j, 'like', k) * k - 0.0012) <0.00000002)
        w0012 = w_prev;
    end
    if(sum(w_prev) <= 0.1 & weGotEm == false)
        weGotEm = true;
        w_final = w_prev;
        t_final = j;
    end
end
t_final = cast(t_final, 'like', k) * k;
%% Graph the nice graph
x = 0:h:L;
x = 1000*x;
plot(x,w0003,x,w0006,x,w0009,x,w0012,x,w_final)
title('Concentration Percent Vs Location at Various Times')
xlabel('Distance (nm)')
ylabel('Concentration (%)')
legend('t = 0.3 ms', 't = 0.6 ms', 't = 0.9 ms', 't = 1.2 ms', 't = 1.5 ms')