(* Created with the Wolfram Language : www.wolfram.com *)
{{2*x4 -> x1 + x2, 2*x3 -> x1 + x4, x1 + x4 -> 2*x1}, 
 SparseArray[Automatic, {5, 5}, 0, {1, {{0, 1, 2, 4, 5, 6}, {{4}, {3}, {2}, 
    {5}, {1}, {3}}}, {kf3, kf1, kr1, kf2, kr3, kr2}}], 
 {{0, 0, 0, 2, 0}, {0, 0, 2, 0, 0}, {1, 0, 0, 1, 0}, {1, 1, 0, 0, 0}, {2, 0, 
  0, 0, 0}}, 
 {{x1 -> ConditionalExpression[(kf2*kr1^2*x3^4)/
      (kf1^2*kr2*((kf2*kr1*x3^2)/(kf1*kr2))^(3/2)), 
     kr3 > 0 && kr2 > 0 && kr1 > 0 && kf3 > 0 && kf2 > 0 && kf1 > 0 && 
      x3 > 0 && x5 > 0], x2 -> ConditionalExpression[
     (kf1*kr3*((kf2*kr1*x3^2)/(kf1*kr2))^(3/2))/(kf3*kr1*x3^2), 
     kr3 > 0 && kr2 > 0 && kr1 > 0 && kf3 > 0 && kf2 > 0 && kf1 > 0 && 
      x3 > 0 && x5 > 0], x4 -> ConditionalExpression[
     Sqrt[(kf2*kr1*x3^2)/(kf1*kr2)], kr3 > 0 && kr2 > 0 && kr1 > 0 && 
      kf3 > 0 && kf2 > 0 && kf1 > 0 && x3 > 0 && x5 > 0]}}}
