E = [0.02, 0.03, 0.04, 0.05, 0.06, 0.08, 0.1, 0.15, 0.2, 0.3, 0.4, 0.5, 0.6, 0.8, 1, 1.5, 2, 3]
sigma_c = [0.147, 0.142, 0.138, 0.134, 0.130, 0.123, 0.117, 0.106, 0.0968, 0.0843, 0.0756, 0.0689, 0.0637, 0.0561,
           0.0503,
           0.0410, 0.0349, 0.0274]
sigma_t = [83.247, 28.642, 13.338, 7.344, 4.52, 2.093, 1.999, 1.906, 0.9398, 0.3733, 0.2166, 0.1512, 0.1175, 0.0846,
           0.0683,
           0.05117, 0.04513, 0.04172]

# density of Pb
p = 11.34


# Energy must be in MeV
def findSigma(e):
    global E
    assert e < E[len(E) - 1]
    if e <= E[0]:
        return sigma_c[0] * p, sigma_t[0] * p
    start = 0
    end = len(E) - 1
    while True:
        mid = int(round((end - start) / 2)) + start
        E_current = E[mid]
        E_next = E[mid + 1]
        if E_current < e <= E_next:
            mid += 1
            break
        elif e == E_current:
            break
        elif e > E_next:
            start = mid
        elif e < E_current:
            end = mid
    return sigma_c[mid] * p, sigma_t[mid] * p