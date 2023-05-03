(* Created with the Wolfram Language : www.wolfram.com *)
{{x3 + x5 -> x1 + x4, 2*x3 -> x2 + x3, x2 + x3 -> 2*x2}, 
 SparseArray[Automatic, {5, 5}, 0, {1, {{0, 1, 2, 4, 5, 6}, {{5}, {3}, {2}, 
    {4}, {3}, {1}}}, {kf1, kf3, kr3, kf2, kr2, kr1}}], 
 {{0, 0, 1, 0, 1}, {0, 0, 2, 0, 0}, {0, 1, 1, 0, 0}, {0, 2, 0, 0, 0}, {1, 0, 
  0, 1, 0}}, 
 {{x1 -> ConditionalExpression[(kr1*((kf3*x2 - kr2*x2)/(2*kr3) + 
        Sqrt[(kf3^2*x2^2 - 2*kf3*kr2*x2^2 + kr2^2*x2^2 + 4*kf2*kr3*x2^2)/
           kr3^2]/2)*x5)/(kf1*x4), kr2 > 0 && kf3 > 0 && kr3 > 0 && 
      kf2 > 0 && kr1 > 0 && kf1 > 0 && x2 > 0 && x4 > 0 && x5 > 0], 
   x3 -> ConditionalExpression[(kf3*x2 - kr2*x2)/(2*kr3) + 
      Sqrt[(kf3^2*x2^2 - 2*kf3*kr2*x2^2 + kr2^2*x2^2 + 4*kf2*kr3*x2^2)/kr3^2]/
       2, kr2 > 0 && kf3 > 0 && kr3 > 0 && kf2 > 0 && kr1 > 0 && kf1 > 0 && 
      x2 > 0 && x4 > 0 && x5 > 0]}}}
