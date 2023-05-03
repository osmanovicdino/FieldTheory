(* Created with the Wolfram Language : www.wolfram.com *)
{{2*x5 -> x1 + x5, x2 + x4 -> 2*x1, x1 + x2 -> 2*x1}, 
 SparseArray[Automatic, {5, 5}, 0, {1, {{0, 1, 2, 3, 4, 6}, {{3}, {5}, {1}, 
    {5}, {2}, {4}}}, {kf2, kf1, kr2, kf3, kr1, kr3}}], 
 {{0, 0, 0, 0, 2}, {0, 1, 0, 1, 0}, {1, 0, 0, 0, 1}, {1, 1, 0, 0, 0}, {2, 0, 
  0, 0, 0}}, {{x1 -> ConditionalExpression[(kf3*kr1*x4)/(kf1*kr3), 
     kr3 > 0 && kr2 > 0 && kr1 > 0 && kf3 > 0 && kf2 > 0 && kf1 > 0 && 
      x3 > 0 && x4 > 0], x2 -> ConditionalExpression[
     ((2*kf3^2*kr1^2*x4^2)/(kf1*kr3^2) + (kf3^3*kr1^2*x4^2)/(kf1^2*kr3^2))/
      (2*kr1*x4 + (kf3*kr1*x4)/kf1), kr3 > 0 && kr2 > 0 && kr1 > 0 && 
      kf3 > 0 && kf2 > 0 && kf1 > 0 && x3 > 0 && x4 > 0], 
   x5 -> ConditionalExpression[(kf2*kf3*kr1*x4)/(kf1*kr2*kr3), 
     kr3 > 0 && kr2 > 0 && kr1 > 0 && kf3 > 0 && kf2 > 0 && kf1 > 0 && 
      x3 > 0 && x4 > 0]}}}
