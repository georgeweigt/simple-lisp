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
(setq A (spacetime-dot a a))
(setq B (sum
  (power at 2)
  (product -1 (power ax 2))
  (product -1 (power ay 2))
  (product -1 (power az 2))
))
(equal A B)
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
(setq A agamma)
(setq B (sum
  (product at gammat)
  (product -1 ax gammax)
  (product -1 ay gammay)
  (product -1 az gammaz)
))
(equal A B)
; note: gammas are square matrices, use "dot" to multiply
; use "spacetime-dot" to multiply spacetime vectors
(setq A (dot agamma bgamma))
(setq B (sum
  (dot -1 bgamma agamma)
  (dot 2 (spacetime-dot a b) I)
))
(equal A B)
(setq A (dot agamma gamma5))
(setq B (dot -1 gamma5 agamma))
(equal A B)
(setq A (dot gammax agamma gammax))
(setq B (sum agamma (dot 2 ax gammax)))
(equal A B)
(setq A (spacetime-dot gamma gamma))
(setq B (dot 4 I))
(equal A B)
(setq A (spacetime-dot gamma (dot agamma gamma)))
(setq B (dot -2 agamma))
(equal A B)
(setq A (spacetime-dot gamma (dot agamma bgamma gamma)))
(setq B (dot 4 (spacetime-dot a b) I))
(equal A B)
(setq A (spacetime-dot gamma (dot agamma bgamma cgamma gamma)))
(setq B (dot -2 cgamma bgamma agamma))
(equal A B)
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
(define expsin (sum
  (product -1/2 i (yexp (product i arg)))
  (product 1/2 i (yexp (product -1 i arg)))
))
(define expcos (sum
  (product 1/2 (yexp (product i arg)))
  (product 1/2 (yexp (product -1 i arg)))
))
(define expsinh (sum
  (product 1/2 (yexp arg))
  (product -1/2 (yexp (product -1 arg)))
))
(define expcosh (sum
  (product 1/2 (yexp arg))
  (product 1/2 (yexp (product -1 arg)))
))
; same as above but for matrices
(define EXP (prog (tmp count)
  (setq tmp 0)
  (setq count order)
loop
  (setq tmp (dot (power count -1) arg (sum I tmp)))
  (setq count (sum count -1))
  (cond ((greaterp count 0) (goto loop)))
  (return (sum I tmp))
))
(define SIN (sum
  (product -1/2 i (EXP (product i arg)))
  (product 1/2 i (EXP (product -1 i arg)))
))
(define COS (sum
  (product 1/2 (EXP (product i arg)))
  (product 1/2 (EXP (product -1 i arg)))
))
(define SINH (sum
  (product 1/2 (EXP arg))
  (product -1/2 (EXP (product -1 arg)))
))
(define COSH (sum
  (product 1/2 (EXP arg))
  (product 1/2 (EXP (product -1 arg)))
))
; for truncating products of power series
(define POWER (cond
  ((greaterp arg2 order) 0)
  (t (list 'power arg1 arg2))
))
(define truncate (eval (subst 'POWER 'power arg)))
(setq A (EXP (dot 1/2 u gammat gammax)))
(setq B (sum
  (product I (expcosh (product 1/2 u)))
  (dot gammat gammax (expsinh (product 1/2 u)))
))
(equal A B)
(setq A (EXP (dot 1/2 theta gammax gammay)))
(setq B (sum
  (product I (expcos (product 1/2 theta)))
  (dot gammax gammay (expsin (product 1/2 theta)))
))
(equal A B)
(setq A (truncate (dot
  (EXP (dot -1/2 u gammat gammaz))
  gammat
  (EXP (dot 1/2 u gammat gammaz))
)))
(setq B (sum
  (product gammat (expcosh u))
  (product gammaz (expsinh u))
))
(equal A B)
(setq A (truncate (dot
  (EXP (dot -1/2 u gammat gammaz))
  gammaz
  (EXP (dot 1/2 u gammat gammaz))
)))
(setq B (sum
  (product gammaz (expcosh u))
  (product gammat (expsinh u))
))
(equal A B)
(setq A (truncate (dot
  (EXP (dot -1/2 u gammat gammaz))
  gammay
  (EXP (dot 1/2 u gammat gammaz))
)))
(equal A gammay)
(setq A (truncate (dot
  (EXP (dot -1/2 u gammat gammaz))
  gammax
  (EXP (dot 1/2 u gammat gammaz))
)))
(equal A gammax)
