(* Created with the Wolfram Language : www.wolfram.com *)
{{x4 + x5 -> 2*x1, x2 + x5 -> 2*x2, 2*x2 -> x1 + x5}, 
 SparseArray[Automatic, {5, 5}, 0, {1, {{0, 1, 2, 4, 5, 6}, {{5}, {3}, {2}, 
    {4}, {3}, {1}}}, {kf1, kf3, kr3, kf2, kr2, kr1}}], 
 {{0, 0, 0, 1, 1}, {0, 1, 0, 0, 1}, {0, 2, 0, 0, 0}, {1, 0, 0, 0, 1}, {2, 0, 
  0, 0, 0}}, {{x1 -> ConditionalExpression[(kr2*kr3*x2)/(kf2*kf3), 
     kr3 > 0 && kr2 > 0 && kr1 > 0 && kf3 > 0 && kf2 > 0 && kf1 > 0 && 
      x2 > 0 && x3 > 0], x4 -> ConditionalExpression[
     (kf1*kr2^2*kr3^3*x2)/(kf2^2*kf3^3*kr1), kr3 > 0 && kr2 > 0 && kr1 > 0 && 
      kf3 > 0 && kf2 > 0 && kf1 > 0 && x2 > 0 && x3 > 0], 
   x5 -> ConditionalExpression[(kf3*x2)/kr3, kr3 > 0 && kr2 > 0 && kr1 > 0 && 
      kf3 > 0 && kf2 > 0 && kf1 > 0 && x2 > 0 && x3 > 0]}}}
