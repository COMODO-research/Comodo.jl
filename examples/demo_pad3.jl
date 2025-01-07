using Comodo

A = rand(50,30,20)

B = pad3(A; padAmount = 2, padValue=1.0)