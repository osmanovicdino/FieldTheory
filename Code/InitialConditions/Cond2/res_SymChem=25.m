(* Created with the Wolfram Language : www.wolfram.com *)
{{x4 + x5 -> 2*x2, 2*x3 -> x2 + x4, 2*x2 -> x1 + x4}, 
 SparseArray[Automatic, {5, 5}, 0, {1, {{0, 1, 2, 3, 5, 6}, {{4}, {3}, {2}, 
    {1}, {5}, {4}}}, {kf3, kf2, kr2, kr3, kf1, kr1}}], 
 {{0, 0, 0, 1, 1}, {0, 0, 2, 0, 0}, {0, 1, 0, 1, 0}, {0, 2, 0, 0, 0}, {1, 0, 
  0, 1, 0}}, {{x1 -> ConditionalExpression[(kr1*kr2^2*x3^4)/(kf1*kf2^2*x4^3), 
     kr3 > 0 && kr2 > 0 && kr1 > 0 && kf3 > 0 && kf2 > 0 && kf1 > 0 && 
      x3 > 0 && x4 > 0], x2 -> ConditionalExpression[(kr2*x3^2)/(kf2*x4), 
     kr3 > 0 && kr2 > 0 && kr1 > 0 && kf3 > 0 && kf2 > 0 && kf1 > 0 && 
      x3 > 0 && x4 > 0], x5 -> ConditionalExpression[
     (kf3*kr2^2*x3^4)/(kf2^2*kr3*x4^3), kr3 > 0 && kr2 > 0 && kr1 > 0 && 
      kf3 > 0 && kf2 > 0 && kf1 > 0 && x3 > 0 && x4 > 0]}}}
