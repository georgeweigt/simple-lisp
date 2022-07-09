; metric tensor

(setq gtt (product -1 (xi r)))
(setq grr (power (xi r) -1))
(setq gthetatheta (power r 2))
(setq gphiphi (product (power r 2) (power (sin theta) 2)))

(setq gdd (sum
  (product gtt (tensor t t))
  (product grr (tensor r r))
  (product gthetatheta (tensor theta theta))
  (product gphiphi (tensor phi phi))
))

(setq g (determinant gdd t r theta phi))
(setq guu (product (power g -1) (adjunct gdd t r theta phi)))

; connection coefficients

(define gradient (sum
  (product (derivative arg t)     (tensor t))
  (product (derivative arg r)     (tensor r))
  (product (derivative arg theta) (tensor theta))
  (product (derivative arg phi)   (tensor phi))
))

(setq gddd (gradient gdd))

(setq GAMDDD (product 1/2 (sum
  gddd
  (transpose 2 3 gddd)
  (product -1 (transpose 1 2 (transpose 2 3 gddd))) ; transpose bc,a to (,a)bc
)))

(setq GAMUDD (dot guu GAMDDD)) ; raise first index

; riemann tensor

(setq GAMUDDD (gradient GAMUDD))

(setq GAMGAM (contract 2 4 (product GAMUDD GAMUDD)))

(setq RUDDD (sum
  (transpose 3 4 GAMUDDD)
  (product -1 GAMUDDD)
  (transpose 2 3 GAMGAM)
  (product -1 (transpose 3 4 (transpose 2 3 GAMGAM)))
))

; ricci tensor

(setq RDD (contract 1 3 RUDDD))

; ricci scalar

(setq R (contract 1 2 (dot guu RDD)))

; einstein tensor

(setq GDD (sum RDD (product -1/2 gdd R)))

"Einstein tensor"
(printcars GDD)

; cancel denominators to simplify

(setq GDD (product GDD r (power r 2) (xi r)))

; replace (xi r) with Schwarzschild metric

(define xi (sum 1 (product -2 M (power r -1))))

(setq GDD (eval GDD))

"Does Einstein tensor vanish for Schwarzschild metric?"
(equal GDD 0)
