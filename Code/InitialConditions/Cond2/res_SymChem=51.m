(* Created with the Wolfram Language : www.wolfram.com *)
{{2*x4 -> x1 + x5, x3 -> 2*x2, x2 + x5 -> 2*x2}, 
 SparseArray[Automatic, {5, 5}, 0, {1, {{0, 1, 2, 3, 5, 6}, {{5}, {4}, {4}, 
    {3}, {2}, {1}}}, {kf3, kf2, kf1, kr1, kr2, kr3}}], 
 {{0, 0, 0, 2, 0}, {0, 0, 1, 0, 0}, {0, 1, 0, 0, 1}, {0, 2, 0, 0, 0}, {1, 0, 
  0, 0, 1}}, {{x1 -> ConditionalExpression[(kr1*kr3*x4^2)/(kf1*kf3*x2), 
     kr3 > 0 && kr2 > 0 && kr1 > 0 && kf3 > 0 && kf2 > 0 && kf1 > 0 && 
      x2 > 0 && x4 > 0], x3 -> ConditionalExpression[(kf2*x2^2)/kr2, 
     kr3 > 0 && kr2 > 0 && kr1 > 0 && kf3 > 0 && kf2 > 0 && kf1 > 0 && 
      x2 > 0 && x4 > 0], x5 -> ConditionalExpression[(kf1*x2)/kr1, 
     kr3 > 0 && kr2 > 0 && kr1 > 0 && kf3 > 0 && kf2 > 0 && kf1 > 0 && 
      x2 > 0 && x4 > 0]}}}
