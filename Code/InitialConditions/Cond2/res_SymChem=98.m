(* Created with the Wolfram Language : www.wolfram.com *)
{{2*x3 -> x1, x2 + x5 -> x1 + x3, x1 + x5 -> x1 + x3}, 
 SparseArray[Automatic, {5, 5}, 0, {1, {{0, 1, 2, 3, 4, 6}, {{3}, {5}, {1}, 
    {5}, {2}, {4}}}, {kf2, kf1, kr2, kf3, kr1, kr3}}], 
 {{0, 0, 2, 0, 0}, {0, 1, 0, 0, 1}, {1, 0, 0, 0, 0}, {1, 0, 0, 0, 1}, {1, 0, 
  1, 0, 0}}, {{x1 -> ConditionalExpression[(2*kr2*x3^2 - (kf1*kr2*x3^3)/kf2)/
      (2*kf2 - kf1*x3), (kf3 > 0 && kf1 > 0 && kr3 > 0 && kr2 > 0 && 
       kr1 > 0 && kf2 > 0 && x3 > (2*kf2)/kf1 && x4 > 0) || 
      (kf3 > 0 && kf1 > 0 && kr3 > 0 && kr2 > 0 && kr1 > 0 && kf2 > 0 && 
       Inequality[0, Less, x3, Less, (2*kf2)/kf1] && x4 > 0)], 
   x2 -> ConditionalExpression[(kf1*kr2*kr3*x3^2)/(kf2*kf3*kr1), 
     (kf3 > 0 && kf1 > 0 && kr3 > 0 && kr2 > 0 && kr1 > 0 && kf2 > 0 && 
       x3 > (2*kf2)/kf1 && x4 > 0) || (kf3 > 0 && kf1 > 0 && kr3 > 0 && 
       kr2 > 0 && kr1 > 0 && kf2 > 0 && Inequality[0, Less, x3, Less, 
        (2*kf2)/kf1] && x4 > 0)], x5 -> ConditionalExpression[(kf3*x3)/kr3, 
     (kf3 > 0 && kf1 > 0 && kr3 > 0 && kr2 > 0 && kr1 > 0 && kf2 > 0 && 
       x3 > (2*kf2)/kf1 && x4 > 0) || (kf3 > 0 && kf1 > 0 && kr3 > 0 && 
       kr2 > 0 && kr1 > 0 && kf2 > 0 && Inequality[0, Less, x3, Less, 
        (2*kf2)/kf1] && x4 > 0)]}, 
  {x1 -> ConditionalExpression[(2*(-2*kf2 - (2*kf2*kf3)/kf1)^2*kr2)/
      ((-kf1 - kf3)^2*((-2*kf2*kf3)/kf1 + (kf1*(-2*kf2 - (2*kf2*kf3)/kf1))/
         (-kf1 - kf3) + (kf3*(-2*kf2 - (2*kf2*kf3)/kf1))/(-kf1 - kf3))), 
     kr3 > 0 && kr2 > 0 && kr1 > 0 && kf3 > 0 && kf2 > 0 && kf1 > 0 && 
      x4 > 0], x2 -> ConditionalExpression[
     (kf1*(-2*kf2 - (2*kf2*kf3)/kf1)^2*kr2*kr3)/(kf2*(-kf1 - kf3)^2*kf3*kr1), 
     kr3 > 0 && kr2 > 0 && kr1 > 0 && kf3 > 0 && kf2 > 0 && kf1 > 0 && 
      x4 > 0], x3 -> ConditionalExpression[(-2*kf2 - (2*kf2*kf3)/kf1)/
      (-kf1 - kf3), kr3 > 0 && kr2 > 0 && kr1 > 0 && kf3 > 0 && kf2 > 0 && 
      kf1 > 0 && x4 > 0], x5 -> ConditionalExpression[(2*kf2*kf3)/(kf1*kr3), 
     kr3 > 0 && kr2 > 0 && kr1 > 0 && kf3 > 0 && kf2 > 0 && kf1 > 0 && 
      x4 > 0]}}}
