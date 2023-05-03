(* Created with the Wolfram Language : www.wolfram.com *)
{{x3 + x4 -> 2*x1, x2 + x5 -> x1 + x3, 2*x2 -> x1 + x3}, 
 SparseArray[Automatic, {5, 5}, 0, {1, {{0, 1, 2, 3, 5, 6}, {{5}, {4}, {4}, 
    {3}, {2}, {1}}}, {kf3, kf2, kf1, kr1, kr2, kr3}}], 
 {{0, 0, 1, 1, 0}, {0, 1, 0, 0, 1}, {0, 2, 0, 0, 0}, {1, 0, 1, 0, 0}, {2, 0, 
  0, 0, 0}}, {{x1 -> ConditionalExpression[(2*kr1*x2^2 + (kf2*kr1*x2^2)/kf1)/
      (2*kf1*x3 + kf2*x3), kf2 > 0 && kf1 > 0 && kr3 > 0 && kr2 > 0 && 
      kr1 > 0 && kf3 > 0 && x2 > 0 && x3 > 0], 
   x4 -> ConditionalExpression[(-(kr1*x2^2) - (kf2*kr1*x2^2)/kf1 + 
       (2*kf3*(2*kr1*x2^2 + (kf2*kr1*x2^2)/kf1)^2)/(2*kf1*x3 + kf2*x3)^2 + 
       (kf1*(2*kr1*x2^2 + (kf2*kr1*x2^2)/kf1)*x3)/(2*kf1*x3 + kf2*x3) + 
       (kf2*(2*kr1*x2^2 + (kf2*kr1*x2^2)/kf1)*x3)/(2*kf1*x3 + kf2*x3))/
      (2*kr3*x3), kf2 > 0 && kf1 > 0 && kr3 > 0 && kr2 > 0 && kr1 > 0 && 
      kf3 > 0 && x2 > 0 && x3 > 0], 
   x5 -> ConditionalExpression[(kf2*kr1*x2)/(kf1*kr2), 
     kf2 > 0 && kf1 > 0 && kr3 > 0 && kr2 > 0 && kr1 > 0 && kf3 > 0 && 
      x2 > 0 && x3 > 0]}}}
