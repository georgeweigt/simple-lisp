; From "Quantum Electrodynamics" by Richard P. Feynman
; pp. 40-43
; generic spacetime vectors a, b and c
(setq a (sum
  (product at (tensor t))
  (product ax (tensor x))
  (product ay (tensor y))
  (product az (tensor z))
))
(setq b (sum
  (product bt (tensor t))
  (product bx (tensor x))
  (product by (tensor y))
  (product bz (tensor z))
))
(setq c (sum
  (product ct (tensor t))
  (product cx (tensor x))
  (product cy (tensor y))
  (product cz (tensor z))
))
; define this function for multiplying spactime vectors
; how it works: (dot arg1 (tensor t)) picks off the t'th element, etc.
; the -1's are for the spacetime metric
(define spacetime-dot (sum
  (dot (dot arg1 (tensor t)) (dot arg2 (tensor t)))
  (dot -1 (dot arg1 (tensor x)) (dot arg2 (tensor x)))
  (dot -1 (dot arg1 (tensor y)) (dot arg2 (tensor y)))
  (dot -1 (dot arg1 (tensor z)) (dot arg2 (tensor z)))
))
(setq tmp1 (spacetime-dot a a))
(setq tmp2 (sum
  (power at 2)
  (product -1 (power ax 2))
  (product -1 (power ay 2))
  (product -1 (power az 2))
))
(equal tmp1 tmp2)
(setq I (sum
  (tensor t t)
  (tensor x x)
  (tensor y y)
  (tensor z z)
))
(setq gammat (sum
  (tensor t t)
  (tensor x x)
  (product -1 (tensor y y))
  (product -1 (tensor z z))
))
(setq gammax (sum
  (tensor t z)
  (tensor x y)
  (product -1 (tensor y x))
  (product -1 (tensor z t))
))
(setq gammay (sum
  (product -1 i (tensor t z))
  (product i (tensor x y))
  (product i (tensor y x))
  (product -1 i (tensor z t))
))
(setq gammaz (sum
  (tensor t y)
  (product -1 (tensor x z))
  (product -1 (tensor y t))
  (tensor z x)
))
(equal (dot gammat gammat) I)
(equal (dot gammax gammax) (dot -1 I))
(equal (dot gammay gammay) (dot -1 I))
(equal (dot gammaz gammaz) (dot -1 I))
(setq gamma5 (dot gammax gammay gammaz gammat))
(equal (dot gamma5 gamma5) (dot -1 I))
; gamma is a "vector" of dirac matrices
(setq gamma (sum
  (product gammat (tensor t))
  (product gammax (tensor x))
  (product gammay (tensor y))
  (product gammaz (tensor z))
))
(setq agamma (spacetime-dot a gamma))
(setq bgamma (spacetime-dot b gamma))
(setq cgamma (spacetime-dot c gamma))
(setq tmp1 agamma)
(setq tmp2 (sum
  (product at gammat)
  (product -1 ax gammax)
  (product -1 ay gammay)
  (product -1 az gammaz)
))
(equal tmp1 tmp2)
; note: gammas are square matrices, use "dot" to multiply
; use "spacetime-dot" to multiply spacetime vectors
(setq tmp1 (dot agamma bgamma))
(setq tmp2 (sum
  (dot -1 bgamma agamma)
  (dot 2 (spacetime-dot a b) I)
))
(equal tmp1 tmp2)
(setq tmp1 (dot agamma gamma5))
(setq tmp2 (dot -1 gamma5 agamma))
(equal tmp1 tmp2)
(setq tmp1 (dot gammax agamma gammax))
(setq tmp2 (sum agamma (dot 2 ax gammax)))
(equal tmp1 tmp2)
(setq tmp1 (spacetime-dot gamma gamma))
(setq tmp2 (dot 4 I))
(equal tmp1 tmp2)
(setq tmp1 (spacetime-dot gamma (dot agamma gamma)))
(setq tmp2 (dot -2 agamma))
(equal tmp1 tmp2)
(setq tmp1 (spacetime-dot gamma (dot agamma bgamma gamma)))
(setq tmp2 (dot 4 (spacetime-dot a b) I))
(equal tmp1 tmp2)
(setq tmp1 (spacetime-dot gamma (dot agamma bgamma cgamma gamma)))
(setq tmp2 (dot -2 cgamma bgamma agamma))
(equal tmp1 tmp2)
; define series approximations for some transcendental functions
; for 32-bit integers, overflow occurs for powers above 5
(define order 5)
(define yexp (prog (tmp count)
  (setq tmp 0)
  (setq count order)
loop
  (setq tmp (product (power count -1) arg (sum 1 tmp)))
  (setq count (sum count -1))
  (cond ((greaterp count 0) (goto loop)))
  (return (sum 1 tmp))
))
(define ysin (sum
  (product -1/2 i (yexp (product i arg)))
  (product 1/2 i (yexp (product -1 i arg)))
))
(define ycos (sum
  (product 1/2 (yexp (product i arg)))
  (product 1/2 (yexp (product -1 i arg)))
))
(define ysinh (sum
  (product 1/2 (yexp arg))
  (product -1/2 (yexp (product -1 arg)))
))
(define ycosh (sum
  (product 1/2 (yexp arg))
  (product 1/2 (yexp (product -1 arg)))
))
; same as above but for matrices
(define YEXP (prog (tmp count)
  (setq tmp 0)
  (setq count order)
loop
  (setq tmp (dot (power count -1) arg (sum I tmp)))
  (setq count (sum count -1))
  (cond ((greaterp count 0) (goto loop)))
  (return (sum I tmp))
))
(define YSIN (sum
  (product -1/2 i (YEXP (product i arg)))
  (product 1/2 i (YEXP (product -1 i arg)))
))
(define YCOS (sum
  (product 1/2 (YEXP (product i arg)))
  (product 1/2 (YEXP (product -1 i arg)))
))
(define YSINH (sum
  (product 1/2 (YEXP arg))
  (product -1/2 (YEXP (product -1 arg)))
))
(define YCOSH (sum
  (product 1/2 (YEXP arg))
  (product 1/2 (YEXP (product -1 arg)))
))
; for truncating products of power series
(define POWER (cond
  ((greaterp arg2 order) 0)
  (t (list 'power arg1 arg2))
))
(define truncate (eval (subst 'POWER 'power arg)))
(setq tmp1 (YEXP (dot 1/2 u gammat gammax)))
(setq tmp2 (sum
  (product I (ycosh (product 1/2 u)))
  (dot gammat gammax (ysinh (product 1/2 u)))
))
(equal tmp1 tmp2)
(setq tmp1 (YEXP (dot 1/2 theta gammax gammay)))
(setq tmp2 (sum
  (product I (ycos (product 1/2 theta)))
  (dot gammax gammay (ysin (product 1/2 theta)))
))
(equal tmp1 tmp2)
(setq tmp1 (truncate (dot
  (YEXP (dot -1/2 u gammat gammaz))
  gammat
  (YEXP (dot 1/2 u gammat gammaz))
)))
(setq tmp2 (sum
  (product gammat (ycosh u))
  (product gammaz (ysinh u))
))
(equal tmp1 tmp2)
(setq tmp1 (truncate (dot
  (YEXP (dot -1/2 u gammat gammaz))
  gammaz
  (YEXP (dot 1/2 u gammat gammaz))
)))
(setq tmp2 (sum
  (product gammaz (ycosh u))
  (product gammat (ysinh u))
))
(equal tmp1 tmp2)
(setq tmp1 (truncate (dot
  (YEXP (dot -1/2 u gammat gammaz))
  gammay
  (YEXP (dot 1/2 u gammat gammaz))
)))
(equal tmp1 gammay)
(setq tmp1 (truncate (dot
  (YEXP (dot -1/2 u gammat gammaz))
  gammax
  (YEXP (dot 1/2 u gammat gammaz))
)))
(equal tmp1 gammax)
