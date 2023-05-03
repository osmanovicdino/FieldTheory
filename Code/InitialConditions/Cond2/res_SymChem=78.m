(* Created with the Wolfram Language : www.wolfram.com *)
{{x5 -> 2*x3, x5 -> x1 + x4, x2 + x3 -> x1 + x5}, 
 SparseArray[Automatic, {5, 5}, 0, {1, {{0, 2, 3, 4, 5, 6}, {{5}, {2}, {1}, 
    {4}, {3}, {1}}}, {kf1, kf2, kr2, kf3, kr3, kr1}}], 
 {{0, 0, 0, 0, 1}, {0, 0, 2, 0, 0}, {0, 1, 1, 0, 0}, {1, 0, 0, 0, 1}, {1, 0, 
  0, 1, 0}}, 
 {{x1 -> ConditionalExpression[(-(kr3*x2*x3) - (kf2*kr1*x3^2)/kr2)/
      (-((kf2*kf3*x3^2)/kr2) - (kf2^2*kf3*kr1*x3^3)/(kr2^2*kr3*x2)), 
     kr3 > 0 && kr2 > 0 && kr1 > 0 && kf3 > 0 && kf2 > 0 && kf1 > 0 && 
      x2 > 0 && x3 > 0], x4 -> ConditionalExpression[
     (kf2^2*kf3*kr1*x3^3)/(kf1*kr2^2*kr3*x2), kr3 > 0 && kr2 > 0 && 
      kr1 > 0 && kf3 > 0 && kf2 > 0 && kf1 > 0 && x2 > 0 && x3 > 0], 
   x5 -> ConditionalExpression[(kf2*x3^2)/kr2, kr3 > 0 && kr2 > 0 && 
      kr1 > 0 && kf3 > 0 && kf2 > 0 && kf1 > 0 && x2 > 0 && x3 > 0]}}}
