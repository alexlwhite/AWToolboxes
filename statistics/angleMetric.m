function M = angleMetric(A,B)


if A<0 && B<0 %switch sign and switch A and B
    B0 = B;
    B = -1*A;
    A = -B0;
end


M = atand(A./B);
%apply a correction so that >45 is positive, <45 is
%negatiion, all the way round the circle
M(B>=0) = M(B>=0)-45;
M(B<0) = M(B<0)+135;

