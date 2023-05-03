(* Created with the Wolfram Language : www.wolfram.com *)
{{x4 -> x1, x2 -> x1 + x3, x1 + x5 -> x1 + x3}, 
 SparseArray[Automatic, {5, 5}, 0, {1, {{0, 1, 2, 3, 4, 6}, {{3}, {5}, {1}, 
    {5}, {2}, {4}}}, {kf2, kf1, kr2, kf3, kr1, kr3}}], 
 {{0, 0, 0, 1, 0}, {0, 1, 0, 0, 0}, {1, 0, 0, 0, 0}, {1, 0, 0, 0, 1}, {1, 0, 
  1, 0, 0}}, 
 {{x1 -> ConditionalExpression[(-(kr1*x2) - (kf2*kr1*x2)/(kf1*x3))/
      (-kf2 - kf1*x3), kr3 > 0 && kr2 > 0 && kr1 > 0 && kf3 > 0 && kf2 > 0 && 
      kf1 > 0 && x3 > 0 && x2 > 0], 
   x4 -> ConditionalExpression[(kf2*kr1*x2)/(kf1*kr2*x3), 
     kr3 > 0 && kr2 > 0 && kr1 > 0 && kf3 > 0 && kf2 > 0 && kf1 > 0 && 
      x3 > 0 && x2 > 0], x5 -> ConditionalExpression[
     ((-kf2 - kf1*x3)*(-(kr1*x2) + (kf1*(-(kr1*x2) - (kf2*kr1*x2)/(kf1*x3))*
          x3)/(-kf2 - kf1*x3) + (kf3*(-(kr1*x2) - (kf2*kr1*x2)/(kf1*x3))*x3)/
         (-kf2 - kf1*x3)))/(kr3*(-(kr1*x2) - (kf2*kr1*x2)/(kf1*x3))), 
     kr3 > 0 && kr2 > 0 && kr1 > 0 && kf3 > 0 && kf2 > 0 && kf1 > 0 && 
      x3 > 0 && x2 > 0]}}}
