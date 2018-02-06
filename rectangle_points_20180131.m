function [x,y] = rectangle_points_20180131(start,stop,y_0,width)

x(1) = start;
y(1) = y_0;

x(2) = stop;
y(2) = y_0;

x(3) = stop;
y(3) = y_0+width;

x(4) = start;
y(4) = y_0+width;

x(5) = start;
y(5) = y_0;