(* Created with the Wolfram Language : www.wolfram.com *)
{{x3 + x5 -> x1 + x2, x1 + x5 -> x1 + x4, x1 + x2 -> 2*x1}, 
 SparseArray[Automatic, {5, 5}, 0, {1, {{0, 1, 2, 3, 5, 6}, {{4}, {3}, {2}, 
    {1}, {5}, {4}}}, {kf3, kf2, kr2, kr3, kf1, kr1}}], 
 {{0, 0, 1, 0, 1}, {1, 0, 0, 0, 1}, {1, 0, 0, 1, 0}, {1, 1, 0, 0, 0}, {2, 0, 
  0, 0, 0}}, {{x1 -> ConditionalExpression[(kr1*x2)/kf1, 
     kr3 > 0 && kr2 > 0 && kr1 > 0 && kf3 > 0 && kf2 > 0 && kf1 > 0 && 
      x2 > 0 && x3 > 0], x4 -> ConditionalExpression[
     (kf3*kr1*kr2*x2^2)/(kf1*kf2*kr3*x3), kr3 > 0 && kr2 > 0 && kr1 > 0 && 
      kf3 > 0 && kf2 > 0 && kf1 > 0 && x2 > 0 && x3 > 0], 
   x5 -> ConditionalExpression[(kf3*kr1*x2^2)/(kf1*kr3*x3), 
     kr3 > 0 && kr2 > 0 && kr1 > 0 && kf3 > 0 && kf2 > 0 && kf1 > 0 && 
      x2 > 0 && x3 > 0]}}}
