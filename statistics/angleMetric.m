function M = angleMetric(A,B)

%If both A and B are negative, switch A to B and change the sign of both
bothNeg = A<0 & B<0;
if any(bothNeg)
    B0 = B;
    B(bothNeg) = -1*A(bothNeg);
    A(bothNeg) = -B0(bothNeg);
end

%compute the angle
M = atand(A./B);

%apply a correction so that >45 is positive, <45 is
%negative, all the way round the circle
M(B>=0) = M(B>=0)-45;
M(B<0)  = M(B<0)+135;

