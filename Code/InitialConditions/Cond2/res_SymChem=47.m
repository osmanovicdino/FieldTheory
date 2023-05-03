(* Created with the Wolfram Language : www.wolfram.com *)
{{x4 -> x1 + x5, x4 + x5 -> x2 + x3, x2 + x5 -> x2 + x3}, 
 SparseArray[Automatic, {5, 5}, 0, {1, {{0, 1, 2, 3, 5, 6}, {{5}, {4}, {4}, 
    {3}, {2}, {1}}}, {kf3, kf2, kf1, kr1, kr2, kr3}}], 
 {{0, 0, 0, 1, 0}, {0, 0, 0, 1, 1}, {0, 1, 0, 0, 1}, {0, 1, 1, 0, 0}, {1, 0, 
  0, 0, 1}}, 
 {{x3 -> ConditionalExpression[(kf2*kr1^2*kr3*x2)/(kf1^2*kf3*kr2*x1), 
     kf1 > 0 && kf2 > 0 && kf3 > 0 && kr1 > 0 && kr2 > 0 && kr3 > 0 && 
      x1 > 0 && x2 > 0], x4 -> ConditionalExpression[
     -1/2*(kr1*x2)/kr2 + Sqrt[(kr1^2*kr3*x2^2 + (4*kf2*kr1^2*kr3*x2^2)/kf1 + 
          (4*kf2^2*kr1^2*kr3*x2^2)/kf1^2)/(kr2^2*kr3)]/2, 
     kf1 > 0 && kf2 > 0 && kf3 > 0 && kr1 > 0 && kr2 > 0 && kr3 > 0 && 
      x1 > 0 && x2 > 0], x5 -> ConditionalExpression[
     (kr3*(-1/2*(kr1*x2)/kr2 + Sqrt[(kr1^2*kr3*x2^2 + (4*kf2*kr1^2*kr3*x2^2)/
             kf1 + (4*kf2^2*kr1^2*kr3*x2^2)/kf1^2)/(kr2^2*kr3)]/2))/(kf3*x1), 
     kf1 > 0 && kf2 > 0 && kf3 > 0 && kr1 > 0 && kr2 > 0 && kr3 > 0 && 
      x1 > 0 && x2 > 0]}}}
