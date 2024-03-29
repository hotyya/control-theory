import sympy as sp
from sympy import symbols, Eq, Function

l = symbols("l")

n = 4
n_c = 3

def form_K(c_h):
    K = list()
    for i in range(n):
        k_temp = list()
        for j in range(i):
            k_temp.append(0)
        k_temp.append(1)
        for j in range(n - i - 1):
            k_temp.append(c_h[j])
        K.append(k_temp)

    return K

def form_gamma(c_h, et_pol_coeffs):
    g = list()
    for i in range(n):
        g_temp = list()
        g_temp.append(c_h[i] - et_pol_coeffs[i])
        g.append(g_temp)

    return g

def stabilization_of_a_fully_controlled_system(P, Q, m):
    S = Q
    P_pow = P
    for i in range(sp.sqrt(P.__len__()) - 1):
        S_temp = sp.Matrix(P_pow * Q)
        S = sp.Matrix.hstack(S, S_temp)
        P_pow = P_pow * P

    print("S = ", S)
    if S.rank() == sp.sqrt(P.__len__()):
        print("There is complete controllability")
    else:
        print("There is no complete controllability")

    c_h = ((P.charpoly()).all_coeffs())[1:]
    print("Coefficients of the characteristic polynomial: ", *c_h)

    K = sp.Matrix(form_K(c_h))
    print("Матрица K: ", K)

    et_pol = sp.simplify(sp.Matrix([[m[0], 0, 0, 0], [0, m[1], 0, 0], [0, 0, m[2], 0], [0, 0, 0, m[3]]]).charpoly())
    et_pol_coeffs = et_pol.all_coeffs()[1:]
    print("коэффициенты эталонного полинома: ", *et_pol_coeffs)

    gamma = sp.Matrix(form_gamma(c_h, et_pol_coeffs))
    print("коэффициенты gamma: ", gamma)

    C = gamma.T * ((S * K).inv())
    print("матрица C: ", C)

    return C

def stabilization_of_a_fully_controlled_system_common(P, Q, m):
    S = Q
    P_pow = P
    for i in range(sp.sqrt(P.__len__()) - 2):
        S_temp = sp.Matrix(P_pow * Q)
        S = sp.Matrix.hstack(S, S_temp)
        P_pow = P_pow * P

    print("S = ", S)
    if S.rank() == sp.sqrt(P.__len__()):
        print("There is complete controllability")
    else:
        print("There is no complete controllability")

    T = Q[:, 0]
    k = [0] * 2
    k[0] += 1
    cnt = 0
    for i in range(n_c - 1):
        cnt = i + 1
        if i != 0:
            t_temp = Q[:, i]
            T_pos = sp.Matrix.hstack(T, t_temp)
            if T_pos.rank() == T.rank() + 1:
                T = T_pos
                k[i] += 1
        if T.rank() == sp.sqrt(T.__len__()):
            break
        for j in range(n_c - 1):
            t_temp = P**(j + 1) * Q[:, i]
            T_pos = sp.Matrix.hstack(T, t_temp)
            if T_pos.rank() == T.rank() + 1:
                T = T_pos
                k[i] += 1
            if T.rank() == sp.sqrt(T.__len__()):
                break

    print("матрица T: ", T)

    P_wave = T.inv() * P * T
    print("матрица P с волной: ", P_wave)

    Q_wave = T.inv() * Q
    print("матрица Q с волной: ", Q_wave)

    K = [] * n_c
    for i in range(3):
        k_temp = [0] * 3
        K.append(k_temp)

    gamma = list()
    for i in range(cnt):
        if k[i] == 1:
            if i == 0:
                c_h = P_wave[0, 0]
            else:
                c_h = P_wave[k[i - 1], k[i - 1]]
                K[k[0]][k[0]] = 1

            if i == 0:
                K[0][0] = 1
            else:
                K[k[0]][k[0]] = 1
        else:
            if i == 0:
                P_temp = P_wave[:k[0], :k[0]]
            else:
                P_temp = P_wave[k[0]:, k[0]:]
            c_h = P_temp[:, :n_c]
            c_h = c_h[:, 1]

            if i == 0:
                for j in range(k[0]):
                    for v in range(k[0]):
                        if v == j:
                            K[j][v] = 1
                        elif v > j:
                            K[j][v] = -c_h[k[0] - (v - j)]
            else:
                for j in range(k[1]):
                    for v in range(k[1]):
                        if v == j:
                            K[k[0] + j][k[0] + v] = 1
                        elif v > j:
                            K[k[0] + j][k[0] + v] = -c_h[k[1] - (v - j)]

        et_pol = list()
        if k[i] != 1:
            for j in range(c_h.__len__()):
                temp = [0] * c_h.__len__()
                temp[j] = m[j]
                et_pol.append(temp)
            et_pol = sp.Matrix(et_pol)
            et_pol_coeffs = (et_pol.charpoly()).all_coeffs()[1:]

            m = m[c_h.__len__():]

            temp_gamma = list()
            if i == 1:
                for j in range(k[0]):
                    temp_gamma.append(0)
            for j in range(c_h.__len__()):
                temp_gamma.append(-c_h[c_h.__len__() - j - 1] - et_pol_coeffs[j])
            if i == 0:
                for j in range(k[1]):
                    temp_gamma.append(0)
        else:
            temp_gamma = list()
            if i == 1:
                for j in range(k[0]):
                    temp_gamma.append(0)
            temp_gamma.append(-c_h + m[0])
            if i == 0:
                for j in range(k[1]):
                    temp_gamma.append(0)
            m = m[1:]

        gamma.append(temp_gamma)

    K = sp.Matrix(K)
    print("матрица K: ", K)

    gamma = sp.Matrix(gamma)
    print("матрица коэффициентов gamma: ", gamma)

    C = gamma * ((T * K).inv())
    print("матрица C: ", C)

    return C