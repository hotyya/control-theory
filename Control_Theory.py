import sympy as sp
import numpy as np
from sympy import symbols, Eq, Function
import stabilization as stbl

t = symbols("t")
h = symbols("h")
f_, g = symbols("f g", cls = Function)

def program_control(P, Q, f, T, x0, x1):
    S = Q
    P_pow = P
    for i in range(sp.sqrt(P.__len__()) - 1):
        S_temp = sp.Matrix(P_pow * Q)
        S = sp.Matrix.hstack(S, (-1)**(i + 1) * S_temp)
        P_pow = P_pow * P

    print("S = ", S)
    if S.rank() == sp.sqrt(P.__len__()):
        print("There is complete controllability")
    else:
        print("There is no complete controllability")

    Y = sp.exp(P * t)
    print("Фундаментальная матрица системы: ", Y)

    B = Y.inv() * Q
    print("Матрица B: ", B)

    A = sp.integrate(B * B.T, (t, 0, T))
    print("Матрица A: ", A)

    n = Y.subs(t, T).inv() * x1 - x0 - sp.integrate(Y.inv() * f, (t, 0, T))
    print("n: ", n)

    if S.rank() == sp.sqrt(P.__len__()):
        C = sp.simplify(A.inv() * n)
    else:
        A_check = sp.Matrix.hstack(A, n)
        A = sp.Matrix(A)
        #A.rang() == A_check.rang()
        if True:
            print("Пара точек x0 x1 управляема")
            C = sp.Matrix([-(10 / (2 - 2 * sp.exp(-2))), 0])
        else:
            return 0

    print("Матрица C: ", C)

    #v = sp.Symbol("v")
    #v = sp.solve(sp.integrate(B * v, (t, 0, T)))
    #print("v: ", v)

    U = sp.simplify(B.T * C)
    print("_______________________________________________________")
    print()
    print("Матрица U: ", U)

    return U

def checking_full_observability(P, R):
    S = R.T
    P_pow = P.T
    for i in range(1, sp.sqrt(P.__len__()) - 2):
        S_temp = P_pow * R.T
        S = sp.Matrix.hstack(S, S_temp)
        P_pow = P_pow * P.T

    print(sp.simplify(S.det()))
    s = ""

    return s

def observation(P, R, f, T):
    S = R.T
    P_pow = P.T
    for i in range(sp.sqrt(P.__len__()) - 1):
        S_temp = P_pow * R.T
        S = sp.Matrix.hstack(S, S_temp)
        P_pow = P_pow * P.T

    print("S = ", S)
    if S.rank() == sp.sqrt(P.__len__()):
        print("There is complete observability")
    else:
        print("There is no complete observability")
        return 0

    Y = sp.simplify(sp.exp(P * t))
    print("Фундаментальная матрица системы: ", Y)

    H = R * Y
    print("Матрица H: ", H)

    D = sp.integrate(H.T * H, (t, 0, T))
    print("Матрица D: ", D)

    g = sp.simplify(sp.Matrix([-t + t**2 - 0.5 * t ** 3]) - H * sp.integrate(Y.inv() * f, (t, 0, t)))
    n = sp.integrate(H.T * g, (t, 0, T))
    print("n: ", n)

    x0 = D.inv() * n
    print("x0: ", x0)

    x_t = sp.simplify(Y * (x0 + sp.integrate(Y.inv() * f, (t, 0, t))))
    print()
    print("x_t: ", x_t)

    return x_t

P = sp.Matrix([[1, 2], [0, 1]])
Q = sp.Matrix([2, 0])
f = sp.Matrix([0, 0])
T = 1
x0 = sp.Matrix([3, 4])
x1 = sp.Matrix([sp.exp(1), 4 * sp.exp(1)])

U = program_control(P, Q, f, T, x0, x1)

print("_________________________________________")
print()
print("№2: ")
print()

gamma = sp.Symbol("gamma")
beta = sp.Symbol("beta")
P = sp.Matrix([[0, -1, 0, -1], [0, 0, 0, 0], [0, 0, -1, -1], [-beta, 0, 0, 0]])
R = sp.Matrix([[beta, 0, 0, 0], [0, 1, gamma, 0]])

print(checking_full_observability(P, R))

print()
print("_________________________________")
print()
print("№3: ")
print()

P = sp.Matrix([[0, 1], [0, 0]])
f = sp.Matrix([1, -1])
R = sp.Matrix([[t, -1]])
T = 1

observation(P, R, f, T)

#P = sp.Matrix([[0, 0, 1, 0], [1, 0, 0, 2], [1, 0, 0, 0], [1, 1, 0, 0]])
#Q = sp.Matrix([0, 0, 1, 0])
#m = sp.Matrix([[-1, -2, -3, -4]])

#print("№ 1: ")
#print()
#stbl.stabilization_of_a_fully_controlled_system(P, Q, m)
#print()
#print("_________________________________________")

#print("№ 2: ")
#print()
#P = sp.Matrix([[1, 0, 0], [1, 1, 1], [-1, 0, 2]])
#Q = sp.Matrix([[0, 1], [1, 1], [1, 0]])
#m = sp.Matrix([[-2, -3, -1]])

#stbl.stabilization_of_a_fully_controlled_system_common(P, Q, m)

#P = sp.Matrix([[-1, 0], [0, 3]])
#Q = sp.Matrix([0, 1])
#f = sp.Matrix([0, 0])
#T = 1
#x0 = sp.Matrix([1, 1])
#x1 = sp.Matrix([0, 0])
#print(program_control(P, Q, f, T, x0, x1))