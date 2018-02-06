function [x,y] = arrow_function_20180131(start,stop,w,dw,l)

x(1) = start;
y(1) = 0;

x(2) = start;
y(2) = w;

x(3) = stop-l;
y(3) = w;

x(4) = stop-l;
y(4) = w+dw;

x(5) = stop;
y(5) = w/2;

x(6) = stop-l;
y(6) = -dw;

x(7) = stop-l;
y(7) = 0;

x(8) = start;
y(8) = 0;
