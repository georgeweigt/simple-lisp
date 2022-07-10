; static spherically symmetric metric

(setq gtt (product -1 (power e (product 2 (Phi r)))))
(setq grr (power e (product 2 (Lambda r))))
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

(setq Gtt (product
  (power e (product 2 (Phi r)))
  (power r -2)
  (derivative (product r (sum 1 (product -1 (power e (product -2 (Lambda r)))))) r)
))

(setq Grr (sum
  (product
    (sum 1 (product -1 (power e (product 2 (Lambda r)))))
    (power r -2)
  )
  (product 2 (derivative (Phi r) r) (power r -1))
))

(setq Gthetatheta (product
  (power r 2)
  (power e (product -2 (Lambda r)))
  (sum
    (derivative (derivative (Phi r) r) r)
    (power (derivative (Phi r) r) 2)
    (product (derivative (Phi r) r) (power r -1))
    (product -1 (derivative (Lambda r) r) (derivative (Phi r) r))
    (product -1 (derivative (Lambda r) r) (power r -1))
  )
))

(setq Gphiphi (product Gthetatheta (power (sin theta) 2)))

(setq G (sum
  (product Gtt (tensor t t))
  (product Grr (tensor r r))
  (product Gthetatheta (tensor theta theta))
  (product Gphiphi (tensor phi phi))
))

"Verify Einstein tensor"
(equal GDD G)
