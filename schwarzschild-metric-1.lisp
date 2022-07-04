; uses (gr) function in init.lisp

; coordinate system

(setq x0 t)
(setq x1 r)
(setq x2 theta)
(setq x3 phi)

(setq gtt (product -1 (xi r)))
(setq grr (power (xi r) -1))
(setq gthetatheta (power r 2))
(setq gphiphi (product (power r 2) (power (sin theta) 2)))

; metric tensor

(setq gdd (sum
  (product gtt (tensor t t))
  (product grr (tensor r r))
  (product gthetatheta (tensor theta theta))
  (product gphiphi (tensor phi phi))
))

(gr) ; compute g, guu, GAMUDD, RUDDD, RDD, R, GDD, GUD and GUU

GDD
