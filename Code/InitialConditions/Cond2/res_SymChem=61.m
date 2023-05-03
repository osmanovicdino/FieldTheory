(* Created with the Wolfram Language : www.wolfram.com *)
{{x3 + x5 -> x1, 2*x3 -> 2*x2, x2 + x3 -> 2*x2}, 
 SparseArray[Automatic, {5, 5}, 0, {1, {{0, 1, 2, 3, 5, 6}, {{5}, {4}, {4}, 
    {3}, {2}, {1}}}, {kf3, kf2, kf1, kr1, kr2, kr3}}], 
 {{0, 0, 1, 0, 1}, {0, 0, 2, 0, 0}, {0, 1, 1, 0, 0}, {0, 2, 0, 0, 0}, {1, 0, 
  0, 0, 0}}, 
 {{x1 -> ConditionalExpression[
     (kr3*(-1/4*(kr1*x2)/kr2 + Sqrt[(kr1^2*x2^2 + 8*kf1*kr2*x2^2 + 
            16*kf2*kr2*x2^2)/kr2^2]/4)*x5)/kf3, kf2 > 0 && kf1 > 0 && 
      kr1 > 0 && kr2 > 0 && kr3 > 0 && kf3 > 0 && x2 > 0 && x4 > 0 && 
      x5 > 0], x3 -> ConditionalExpression[-1/4*(kr1*x2)/kr2 + 
      Sqrt[(kr1^2*x2^2 + 8*kf1*kr2*x2^2 + 16*kf2*kr2*x2^2)/kr2^2]/4, 
     kf2 > 0 && kf1 > 0 && kr1 > 0 && kr2 > 0 && kr3 > 0 && kf3 > 0 && 
      x2 > 0 && x4 > 0 && x5 > 0]}}}
