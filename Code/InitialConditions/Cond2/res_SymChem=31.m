(* Created with the Wolfram Language : www.wolfram.com *)
{{x3 + x5 -> x1 + x3, x2 + x4 -> x2 + x3, x2 + x3 -> x1 + x2}, 
 SparseArray[Automatic, {5, 5}, 0, {1, {{0, 1, 2, 4, 5, 6}, {{4}, {3}, {2}, 
    {5}, {1}, {3}}}, {kf3, kf1, kr1, kf2, kr3, kr2}}], 
 {{0, 0, 1, 0, 1}, {0, 1, 0, 1, 0}, {0, 1, 1, 0, 0}, {1, 0, 1, 0, 0}, {1, 1, 
  0, 0, 0}}, 
 {{x1 -> ConditionalExpression[(-((kr1*kr2*x2*x4)/kf1) - 
       (kf3*kr1^2*kr2*x4^2)/(kf1^2*kf2))/(-(kf2*x2) - (kf3*kr1*x4)/kf1), 
     kr3 > 0 && kr2 > 0 && kr1 > 0 && kf3 > 0 && kf2 > 0 && kf1 > 0 && 
      x2 > 0 && x4 > 0], x3 -> ConditionalExpression[(kr1*x4)/kf1, 
     kr3 > 0 && kr2 > 0 && kr1 > 0 && kf3 > 0 && kf2 > 0 && kf1 > 0 && 
      x2 > 0 && x4 > 0], x5 -> ConditionalExpression[
     (kf3*kr1*kr2*x4)/(kf1*kf2*kr3), kr3 > 0 && kr2 > 0 && kr1 > 0 && 
      kf3 > 0 && kf2 > 0 && kf1 > 0 && x2 > 0 && x4 > 0]}}}
