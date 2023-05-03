(* Created with the Wolfram Language : www.wolfram.com *)
{{x5 -> 2*x2, x3 + x4 -> x2 + x5, 2*x2 -> x1 + x3}, 
 SparseArray[Automatic, {5, 5}, 0, {1, {{0, 1, 2, 3, 5, 6}, {{4}, {3}, {2}, 
    {1}, {5}, {4}}}, {kf3, kf2, kr2, kr3, kf1, kr1}}], 
 {{0, 0, 0, 0, 1}, {0, 0, 1, 1, 0}, {0, 1, 0, 0, 1}, {0, 2, 0, 0, 0}, {1, 0, 
  1, 0, 0}}, {{x1 -> ConditionalExpression[(kr1*kr3*x5)/(kf1*kf3*x3), 
     kr3 > 0 && kr2 > 0 && kr1 > 0 && kf3 > 0 && kf2 > 0 && kf1 > 0 && 
      x3 > 0 && x5 > 0], x2 -> ConditionalExpression[
     (kr2*x3*Sqrt[(kf2^2*kr3*x5^3)/(kf3*kr2^2*x3^2)])/(kf2*x5), 
     kr3 > 0 && kr2 > 0 && kr1 > 0 && kf3 > 0 && kf2 > 0 && kf1 > 0 && 
      x3 > 0 && x5 > 0], x4 -> ConditionalExpression[
     Sqrt[(kf2^2*kr3*x5^3)/(kf3*kr2^2*x3^2)], kr3 > 0 && kr2 > 0 && 
      kr1 > 0 && kf3 > 0 && kf2 > 0 && kf1 > 0 && x3 > 0 && x5 > 0]}}}
