(* Created with the Wolfram Language : www.wolfram.com *)
{{2*x5 -> 2*x1, x4 + x5 -> x1 + x5, x1 + x5 -> x1 + x2}, 
 SparseArray[Automatic, {5, 5}, 0, {1, {{0, 1, 2, 4, 5, 6}, {{5}, {3}, {2}, 
    {4}, {3}, {1}}}, {kf1, kf3, kr3, kf2, kr2, kr1}}], 
 {{0, 0, 0, 0, 2}, {0, 0, 0, 1, 1}, {1, 0, 0, 0, 1}, {1, 1, 0, 0, 0}, {2, 0, 
  0, 0, 0}}, {{x1 -> ConditionalExpression[(kr3*x4)/kf3, 
     kr3 > 0 && kr2 > 0 && kr1 > 0 && kf3 > 0 && kf2 > 0 && kf1 > 0 && 
      x3 > 0 && x4 > 0], x2 -> ConditionalExpression[
     (kr2*Sqrt[(kf1*kr3^2*x4^2)/(kf3^2*kr1)])/kf2, kr3 > 0 && kr2 > 0 && 
      kr1 > 0 && kf3 > 0 && kf2 > 0 && kf1 > 0 && x3 > 0 && x4 > 0], 
   x5 -> ConditionalExpression[Sqrt[(kf1*kr3^2*x4^2)/(kf3^2*kr1)], 
     kr3 > 0 && kr2 > 0 && kr1 > 0 && kf3 > 0 && kf2 > 0 && kf1 > 0 && 
      x3 > 0 && x4 > 0]}}}
