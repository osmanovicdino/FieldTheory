(* Created with the Wolfram Language : www.wolfram.com *)
{{x5 -> 2*x4, x5 -> 2*x1, x2 -> x1 + x4}, SparseArray[Automatic, {5, 5}, 0, 
  {1, {{0, 2, 3, 4, 5, 6}, {{5}, {2}, {1}, {4}, {3}, {1}}}, 
   {kf1, kf2, kr2, kf3, kr3, kr1}}], {{0, 0, 0, 0, 1}, {0, 0, 0, 2, 0}, {0, 
  1, 0, 0, 0}, {1, 0, 0, 1, 0}, {2, 0, 0, 0, 0}}, 
 {{x1 -> ConditionalExpression[
     (kr3*Sqrt[(kf2*kf3^2*kr1*x4^4)/(kf1*kr2*kr3^2)])/(kf3*x4), 
     kr3 > 0 && kr2 > 0 && kr1 > 0 && kf3 > 0 && kf2 > 0 && kf1 > 0 && 
      x3 > 0 && x4 > 0], x2 -> ConditionalExpression[
     Sqrt[(kf2*kf3^2*kr1*x4^4)/(kf1*kr2*kr3^2)], kr3 > 0 && kr2 > 0 && 
      kr1 > 0 && kf3 > 0 && kf2 > 0 && kf1 > 0 && x3 > 0 && x4 > 0], 
   x5 -> ConditionalExpression[(kf2*x4^2)/kr2, kr3 > 0 && kr2 > 0 && 
      kr1 > 0 && kf3 > 0 && kf2 > 0 && kf1 > 0 && x3 > 0 && x4 > 0]}}}
