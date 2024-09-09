function A=get_element_area(xe,b)

%tensor basis
e1=b(1);
e2=b(2);

% corner coordinates
x1=dot(xe(1),e1);
y1=dot(xe(1),e2);

x2=dot(xe(2),e1);
y2=dot(xe(2),e2);

x3=dot(xe(3),e1);
y3=dot(xe(3),e2);

x4=dot(xe(4),e1);
y4=dot(xe(4),e2);


A=1/2* ( (x1-x3)*(y2-y4) - (x2-x4)*(y1-y3) );

end