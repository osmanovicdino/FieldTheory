(* Created with the Wolfram Language : www.wolfram.com *)
{{x5 -> x4, x5 -> x1 + x2, x2 -> x1}, SparseArray[Automatic, {5, 5}, 0, 
  {1, {{0, 2, 3, 4, 5, 6}, {{5}, {2}, {1}, {4}, {3}, {1}}}, 
   {kf1, kf2, kr2, kf3, kr3, kr1}}], {{0, 0, 0, 0, 1}, {0, 0, 0, 1, 0}, {0, 
  1, 0, 0, 0}, {1, 0, 0, 0, 0}, {1, 1, 0, 0, 0}}, 
 {{x1 -> ConditionalExpression[(-(kr3*x2) - (kf1*kr3*x2^2)/kf3)/
      (-kf3 - kf1*x2), kr3 > 0 && kr2 > 0 && kr1 > 0 && kf3 > 0 && kf2 > 0 && 
      kf1 > 0 && x2 > 0 && x3 > 0], 
   x4 -> ConditionalExpression[((kf1*kr3*x2^2)/kf3 + (kf1*kr2*kr3*x2^2)/
        (kf3*kr1) - (kf1*x2*(-(kr3*x2) - (kf1*kr3*x2^2)/kf3))/
        (-kf3 - kf1*x2))/kf2, kr3 > 0 && kr2 > 0 && kr1 > 0 && kf3 > 0 && 
      kf2 > 0 && kf1 > 0 && x2 > 0 && x3 > 0], 
   x5 -> ConditionalExpression[(kf1*kr3*x2^2)/(kf3*kr1), 
     kr3 > 0 && kr2 > 0 && kr1 > 0 && kf3 > 0 && kf2 > 0 && kf1 > 0 && 
      x2 > 0 && x3 > 0]}}}
