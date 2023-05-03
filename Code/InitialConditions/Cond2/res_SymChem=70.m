(* Created with the Wolfram Language : www.wolfram.com *)
{{x3 + x5 -> x1 + x4, x2 + x4 -> 2*x2, 2*x2 -> x1}, 
 SparseArray[Automatic, {5, 5}, 0, {1, {{0, 1, 2, 4, 5, 6}, {{5}, {3}, {2}, 
    {4}, {3}, {1}}}, {kf1, kf3, kr3, kf2, kr2, kr1}}], 
 {{0, 0, 1, 0, 1}, {0, 1, 0, 1, 0}, {0, 2, 0, 0, 0}, {1, 0, 0, 0, 0}, {1, 0, 
  0, 1, 0}}, 
 {{x1 -> ConditionalExpression[(-((kr2*kr3^2*x4^2)/kf3^2) - 
       (kf1*kr2*kr3^2*x4^3)/(kf2*kf3^2))/(-kf2 - kf1*x4), 
     kr2 > 0 && kf3 > 0 && kr3 > 0 && kr1 > 0 && kf2 > 0 && kf1 > 0 && 
      x4 > 0 && x3 > 0], x2 -> ConditionalExpression[(kr3*x4)/kf3, 
     kr2 > 0 && kf3 > 0 && kr3 > 0 && kr1 > 0 && kf2 > 0 && kf1 > 0 && 
      x4 > 0 && x3 > 0], x5 -> ConditionalExpression[
     (kf1*kr2*kr3^2*x4^3)/(kf2*kf3^2*kr1*x3), kr2 > 0 && kf3 > 0 && 
      kr3 > 0 && kr1 > 0 && kf2 > 0 && kf1 > 0 && x4 > 0 && x3 > 0]}}}
