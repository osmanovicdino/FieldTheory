(* Created with the Wolfram Language : www.wolfram.com *)
{{2*x5 -> 2*x1, x3 -> x1 + x5, 2*x2 -> x1 + x5}, 
 SparseArray[Automatic, {5, 5}, 0, {1, {{0, 1, 2, 3, 5, 6}, {{5}, {4}, {4}, 
    {3}, {2}, {1}}}, {kf3, kf2, kf1, kr1, kr2, kr3}}], 
 {{0, 0, 0, 0, 2}, {0, 0, 1, 0, 0}, {0, 2, 0, 0, 0}, {1, 0, 0, 0, 1}, {2, 0, 
  0, 0, 0}}, 
 {{x1 -> ConditionalExpression[(kr1*x2^2)/
      (kf1*Root[-(kf3*kr1^2*x2^4) + kf1^2*kr3*#1^4 & , 2]), 
     kr3 > 0 && kr2 > 0 && kr1 > 0 && kf3 > 0 && kf2 > 0 && kf1 > 0 && 
      x2 > 0 && x4 > 0], x3 -> ConditionalExpression[
     ((kf2*kr1*x2^2)/kf1 + (2*kf3*kr1^2*x2^4)/
        (kf1^2*Root[-(kf3*kr1^2*x2^4) + kf1^2*kr3*#1^4 & , 2]^2) - 
       2*kr3*Root[-(kf3*kr1^2*x2^4) + kf1^2*kr3*#1^4 & , 2]^2)/kr2, 
     kr3 > 0 && kr2 > 0 && kr1 > 0 && kf3 > 0 && kf2 > 0 && kf1 > 0 && 
      x2 > 0 && x4 > 0], x5 -> ConditionalExpression[
     Root[-(kf3*kr1^2*x2^4) + kf1^2*kr3*#1^4 & , 2], 
     kr3 > 0 && kr2 > 0 && kr1 > 0 && kf3 > 0 && kf2 > 0 && kf1 > 0 && 
      x2 > 0 && x4 > 0]}}}
