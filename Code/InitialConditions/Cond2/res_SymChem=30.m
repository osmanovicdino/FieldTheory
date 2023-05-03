(* Created with the Wolfram Language : www.wolfram.com *)
{{x4 -> x1 + x5, x4 + x5 -> x3 + x4, x1 + x5 -> x1 + x3}, 
 SparseArray[Automatic, {5, 5}, 0, {1, {{0, 1, 2, 3, 5, 6}, {{4}, {3}, {2}, 
    {1}, {5}, {4}}}, {kf3, kf2, kr2, kr3, kf1, kr1}}], 
 {{0, 0, 0, 1, 0}, {0, 0, 0, 1, 1}, {0, 0, 1, 1, 0}, {1, 0, 0, 0, 1}, {1, 0, 
  1, 0, 0}}, 
 {{x1 -> ConditionalExpression[
     (-(kr2*x4*x5) + (kf2*x4*(kr1*kr3*x5 + kf3*kr2*x5^2))/
        (kf1*kr3 + kf2*kf3*x5))/(kr1*x5 - (kf1*(kr1*kr3*x5 + kf3*kr2*x5^2))/
        (kf1*kr3 + kf2*kf3*x5)), kr2 > 0 && kr1 > 0 && kf2 > 0 && kf1 > 0 && 
      kr3 > 0 && kf3 > 0 && x2 > 0 && x5 > 0 && x4 > 0], 
   x3 -> ConditionalExpression[(kr1*kr3*x5 + kf3*kr2*x5^2)/
      (kf1*kr3 + kf2*kf3*x5), kr2 > 0 && kr1 > 0 && kf2 > 0 && kf1 > 0 && 
      kr3 > 0 && kf3 > 0 && x2 > 0 && x5 > 0 && x4 > 0]}}}
