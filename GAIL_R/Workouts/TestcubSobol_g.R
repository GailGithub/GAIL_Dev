# Some tests of cubSobol_g.R
answer = cubSobol_g(f = function(x) {sin(x)}) 
answer 
error = answer[1]-(1-cos(1)) 
error

answer = cubSobol_g(f = function(x) {sin(x)},reltol = 0)
answer
error = answer[1]-(1-cos(1))
error

answer = cubSobol_g(f = function(x) {sin(x)},abstol = 1e-8, reltol = 0) 
answer
error = answer[1]-(1-cos(1))
error
warnings()