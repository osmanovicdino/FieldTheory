(* Created with the Wolfram Language : www.wolfram.com *)
{{x4 -> x3, x4 + x5 -> 2*x1, 2*x2 -> 2*x1}, SparseArray[Automatic, {5, 5}, 0, 
  {1, {{0, 1, 2, 3, 4, 6}, {{3}, {5}, {1}, {5}, {2}, {4}}}, 
   {kf2, kf1, kr2, kf3, kr1, kr3}}], {{0, 0, 0, 1, 0}, {0, 0, 0, 1, 1}, {0, 
  0, 1, 0, 0}, {0, 2, 0, 0, 0}, {2, 0, 0, 0, 0}}, 
 {{x2 -> ConditionalExpression[Sqrt[(kf3*x1^2)/kr3], 
     kr3 > 0 && kr2 > 0 && kr1 > 0 && kf3 > 0 && kf2 > 0 && kf1 > 0 && 
      x1 > 0 && x5 > 0], x3 -> ConditionalExpression[
     (kf1*kr2*x1^2)/(kf2*kr1*x5), kr3 > 0 && kr2 > 0 && kr1 > 0 && kf3 > 0 && 
      kf2 > 0 && kf1 > 0 && x1 > 0 && x5 > 0], 
   x4 -> ConditionalExpression[(kf1*x1^2)/(kr1*x5), 
     kr3 > 0 && kr2 > 0 && kr1 > 0 && kf3 > 0 && kf2 > 0 && kf1 > 0 && 
      x1 > 0 && x5 > 0]}}}
