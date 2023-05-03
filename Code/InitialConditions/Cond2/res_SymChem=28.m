(* Created with the Wolfram Language : www.wolfram.com *)
{{x4 + x5 -> x3 + x4, x3 + x5 -> x1, 2*x3 -> x1}, 
 SparseArray[Automatic, {5, 5}, 0, {1, {{0, 1, 2, 3, 4, 6}, {{3}, {5}, {1}, 
    {5}, {2}, {4}}}, {kf2, kf1, kr2, kf3, kr1, kr3}}], 
 {{0, 0, 0, 1, 1}, {0, 0, 1, 0, 1}, {0, 0, 1, 1, 0}, {0, 0, 2, 0, 0}, {1, 0, 
  0, 0, 0}}, 
 {{x1 -> ConditionalExpression[(-(kr3*x3^2) - kr1*x3*x5)/(-kf1 - kf3), 
     (kr3 > 0 && kr2 > 0 && kr1 > 0 && kf3 > 0 && kf2 > 0 && 
       kf1 > (kf2*kf3*kr1)/(kr2*kr3) && x2 > 0 && x3 > 0 && 
       Inequality[(kf2*x3)/kr2, Less, x5, Less, (kf1*kr3*x3)/(kf3*kr1)]) || 
      (kr3 > 0 && kr2 > 0 && kr1 > 0 && kf3 > 0 && kf2 > 0 && 
       Inequality[0, Less, kf1, Less, (kf2*kf3*kr1)/(kr2*kr3)] && x2 > 0 && 
       x3 > 0 && Inequality[(kf1*kr3*x3)/(kf3*kr1), Less, x5, Less, 
        (kf2*x3)/kr2])], x4 -> ConditionalExpression[
     (2*kr3*x3^2 + kr1*x3*x5 - (kf1*(-(kr3*x3^2) - kr1*x3*x5))/(-kf1 - kf3) - 
       (2*kf3*(-(kr3*x3^2) - kr1*x3*x5))/(-kf1 - kf3))/(-(kf2*x3) + kr2*x5), 
     (kr3 > 0 && kr2 > 0 && kr1 > 0 && kf3 > 0 && kf2 > 0 && 
       kf1 > (kf2*kf3*kr1)/(kr2*kr3) && x2 > 0 && x3 > 0 && 
       Inequality[(kf2*x3)/kr2, Less, x5, Less, (kf1*kr3*x3)/(kf3*kr1)]) || 
      (kr3 > 0 && kr2 > 0 && kr1 > 0 && kf3 > 0 && kf2 > 0 && 
       Inequality[0, Less, kf1, Less, (kf2*kf3*kr1)/(kr2*kr3)] && x2 > 0 && 
       x3 > 0 && Inequality[(kf1*kr3*x3)/(kf3*kr1), Less, x5, Less, 
        (kf2*x3)/kr2])]}}}
