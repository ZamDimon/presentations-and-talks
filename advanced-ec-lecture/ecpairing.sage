import json

# --- Below we specify the parameters of the curve --- #

# Defining the base prime field
q = Integer(21888242871839275222246405745257275088696311157297823662689037894645226208583) # EC group order
Fq = GF(q) 

# r is taken from https://hackmd.io/@jpw/bn254
k = Integer(12) # Embedding degree
t = Integer(4965661367192848881)
r = Integer(21888242871839275222246405745257275088548364400416034343698204186575808495617)
e = (q^(12)-1)/r

# Making sure parameters are correctly defined
# See https://eprint.iacr.org/2010/354.pdf, Equation 1 for details.
assert q == 36*t**4 + 36*t**3 + 24*t**2 + 6*t + 1
assert r == 36*t**4 + 36*t**3 + 18*t**2 + 6*t + 1

# Defining the extensions
# Fq2...
K2.<x> = PolynomialRing(Fq)
Fq2.<u> = Fq.extension(x^2+1)

# Fq6...
K6.<y> = PolynomialRing(Fq2)
Fq6.<v> = Fq2.extension(y^3 - (u+9))

# Defining the Fq12 is a bit more tricky...
p = Fq.characteristic()
Fq12.<G> = GF(p^12)

i = sqrt(Fq12(-1))
R12.<Y> = PolynomialRing(Fq12)

j = (Y^3 - (i+9)).roots(multiplicities=False)[0]
w = sqrt(j)

P = w.minpoly()
Fq12.<W> = GF(p^12, modulus=P)

def c0c3c4_to_fq12(c0: Fq2, c3: Fq2, c4: Fq2) -> Fq12:
    return c0[0] + c0[1]*(W^6-9) + (c3[0]+c3[1]*(W^6-9))*W + (c4[0]+c4[1]*(W^6-9))*W^3

# Defining the G1 Curve and its generator
G1 = EllipticCurve(Fq, [0, 3])
G1_GEN = G1(1, 2)

# Defining the G2 Curve
b = 3 / (u + 9)
G2 = EllipticCurve(Fq2, [0, b])
G2_GEN = G2(10857046999023057135944570762232829481370756359578518086990519993285655852781+
            11559732032986387107991004021392285783925812861821192530917403151452391805634*u,
            8495653923123431417604973247489272438418190587263600148770280649306958101930+
            4082367875863433681332203403145435568316851327593401208105741076214120093531*u)

# Converts a tuple (X : Y : Z) from Fq2^3 to a point in G2 
# using Jacobian coordinates
def tuple_to_g2(t: tuple[Fq2, Fq2, Fq2]) -> G2:
    return G2(t[0]/t[2]^2, t[1]/t[2]^3)

# Helper debugging functions
g1_point_to_dictionary = lambda point : {
    'x': str(point[0]),
    'y': str(point[1])
}
g2_point_to_dictionary = lambda point : {
    'x': {
        'c0': str(point[0][0]), 
        'c1': str(point[0][1])
    }, 
    'y': {
        'c0': str(point[1][0]), 
        'c1': str(point[1][1])
    }
}

# Some coefficients for easier life
SIX_U_PLUS_TWO_WNAF = [
    0, 0, 0, 1, 0, 1, 0, -1, 
    0, 0, 1, -1, 0, 0, 1, 0, 
    0, 1, 1, 0, -1, 0, 0, 1, 
    0, -1, 0, 0, 0, 0, 1, 1, 
    1, 0, 0, -1, 0, 0, 1, 0, 
    0, 0, 0, 0, -1, 0, 0, 1, 
    1, 0, 0, -1, 0, 0, 0, 1, 
    1, 0, -1, 0, 0, 1, 0, 1, 1
]

# Converts the Montomery form represented by 4 64-bit limbs to an integer in Fq
def from_libms(limbs):
    montomery = limbs[0] | (limbs[1] << 64) | (limbs[2] << 128) | (limbs[3] << 192)
    return Fq(montomery) * Fq(2^(-256))

# This is for the last step of Miller loop
FROBENIUS_COEFF_FQ6_C1_1 = from_libms([
    0xb5773b104563ab30,
    0x347f91c8a9aa6454,
    0x7a007127242e0991,
    0x1956bcd8118214ec,
]) + from_libms([
    0x6e849f1ea0aa4757, 
    0xaa1c7b6d89f89141, 
    0xb6e713cdfae0ca3a, 
    0x26694fbb4e82ebc3,
])*u
assert FROBENIUS_COEFF_FQ6_C1_1 == (9+u)**((q-1)/3), 'FROBENIUS_COEFF_FQ6_C1_1 is not correct!'

# (9+u)**((q-1)/2)
XI_TO_Q_MINUS_1_OVER_2 = from_libms([
    0xe4bbdd0c2936b629, 
    0xbb30f162e133bacb, 
    0x31a9d1b6f9645366, 
    0x253570bea500f8dd,
]) + from_libms([
    0xa1d77ce45ffe77c7, 
    0x07affd117826d1db, 
    0x6d16bd27bb7edc6b, 
    0x2c87200285defecc,
])*u
assert XI_TO_Q_MINUS_1_OVER_2 == (9+u)**((q-1)/2), 'Non-XI_TO_Q_MINUS_1_OVER_2 is not correct!'

# (9+u)**((q^2-1)/3)
FROBENIUS_COEFF_FQ6_C1_2 = from_libms([
    0x3350c88e13e80b9c,
    0x7dce557cdb5e56b9,
    0x6001b4b8b615564a,
    0x2682e617020217e0,
]) + from_libms([
    0x0, 
    0x0, 
    0x0, 
    0x0,
])*u
assert FROBENIUS_COEFF_FQ6_C1_2 == (9+u)**((q^2-1)/3), 'FROBENIUS_COEFF_FQ6_C1_2 is not correct!'

# --- Line functions tested ---
# Original implementation from https://eprint.iacr.org/2010/354.pdf

def doubling_step(Q: G2, P: G2):
    X_Q, Y_Q, Z_Q = copy(Q[0]), copy(Q[1]), copy(Q[2])
    x_P, y_P = copy(P[0]), copy(P[1])

    tmp0 = X_Q**2
    tmp1 = Y_Q**2
    tmp2 = tmp1^2
    tmp3 = (tmp1 + X_Q)^2 - tmp0 - tmp2
    tmp3 = 2*tmp3
    tmp4 = 3*tmp0
    tmp6 = X_Q + tmp4
    tmp5 = tmp4^2
    X_T = tmp5 - 2*tmp3
    Z_T = (Y_Q + Z_Q)^2 - tmp1 - Z_Q^2
    Y_T = (tmp3 - X_T) * tmp4 - 8*tmp2
    tmp3 = -2*tmp4*Z_Q^2
    tmp3 = tmp3*x_P
    tmp6 = tmp6^2 - tmp0 - tmp5 - 4*tmp1
    tmp0 = 2*Z_T*Z_Q^2
    tmp0 = tmp0 * y_P

    return (tmp0, tmp3, tmp6), (X_T, Y_T, Z_T)

def addition_step(Q: G2, R: G2, P: G1):
    X_Q, Y_Q, Z_Q = copy(Q[0]), copy(Q[1]), copy(Q[2])
    X_R, Y_R, Z_R = copy(R[0]), copy(R[1]), copy(R[2])
    x_P, y_P = copy(P[0]), copy(P[1])

    t0 = X_Q * Z_R^2
    t1 = (Y_Q + Z_R)^2 - Y_Q^2 - Z_R^2
    t1 = t1 * Z_R^2
    t2 = t0 - X_R
    t3 = t2^2 
    t4 = 4*t3
    t5 = t4 * t2
    t6 = t1 - 2*Y_R
    t9 = t6 * X_Q
    t7 = X_R*t4
    X_T = t6^2 - t5 - 2*t7
    Z_T = (Z_R + t2)^2 - Z_R^2 - t3
    t10 = Y_Q + Z_T
    t8 = (t7 - X_T)*t6
    t0 = 2*Y_R*t5
    Y_T = t8 - t0
    t10 = t10^2 - Y_Q^2 - Z_T^2
    t9 = 2*t9 - t10
    t10 = 2*Z_T*y_P
    t6 = -t6
    t1 = 2*t6*x_P

    return (t10, t1, t9), (X_T, Y_T, Z_T)

def miller_loop(P: G1, Q: G2):
    # --- Gathering coefficients step ---
    T = copy(Q)
    Q_negative = -copy(Q)
    f = Fq12.one()
    for i in reversed(range(1, len(SIX_U_PLUS_TWO_WNAF))):
        if i != len(SIX_U_PLUS_TWO_WNAF) - 1:
            f = f*f

        (c0, c3, c4), T2 = doubling_step(T, P)
        assert tuple_to_g2(T2) == 2*tuple_to_g2(T), 'Doubling step is wrong!'
        f = f * c0c3c4_to_fq12(c0, c3, c4)
        T = T2

        x = SIX_U_PLUS_TWO_WNAF[i-1]
        if x == 1:
            (c0, c3, c4), TQ = addition_step(Q, T, P)
            assert tuple_to_g2(TQ) == tuple_to_g2(T) + tuple_to_g2(Q), 'Addition step is wrong!'
            f = f * c0c3c4_to_fq12(c0, c3, c4)
            T = TQ
        elif x == -1:
            (c0, c3, c4), TQ = addition_step(Q_negative, T, P)
            assert tuple_to_g2(TQ) == tuple_to_g2(T) + Q_negative, 'Addition step is wrong!'
            f = f * c0c3c4_to_fq12(c0, c3, c4)
            T = TQ

    # Some additional steps to finalize the Miller loop...
    # Q1 <- pi_p(Q)
    Q1 = [Q[0], Q[1], Q[2]]
    Q1[0] = Q1[0].conjugate() * FROBENIUS_COEFF_FQ6_C1_1
    Q1[1] = Q1[1].conjugate() * XI_TO_Q_MINUS_1_OVER_2

    # Q2 <- -pi_{p^2}(Q)
    Q2 = [Q[0], Q[1], Q[2]]
    Q2[0] = Q2[0] * FROBENIUS_COEFF_FQ6_C1_2

    # Line evaluation at Q1
    (c0, c3, c4), TQ1 = addition_step(Q1, T, P)
    assert tuple_to_g2(TQ1) == tuple_to_g2(T) + tuple_to_g2(Q1), 'Addition step is wrong!'
    f = f * c0c3c4_to_fq12(c0, c3, c4)
    T = TQ1

    # Line evaluation at Q2
    (c0, c3, c4), TQ2 = addition_step(Q2, T, P)
    assert tuple_to_g2(TQ2) == tuple_to_g2(T) + tuple_to_g2(Q2), 'Addition step is wrong!'
    f = f * c0c3c4_to_fq12(c0, c3, c4)
    
    return f

def pairing(P, Q):
    f = miller_loop(P, Q)
    return f^e # This part is actually the final exponentiation which is usually done differently

# --- Actual demo below --- #

# Defining random elements and points
a = Fq.random_element()
A = a * G1_GEN

b = Fq.random_element()
B = b * G2_GEN

# We are debugging points chosen:
print(f'\n\nPoint on G1: {g1_point_to_dictionary(A)}')
print(f'\n\nPoint on G2: {g2_point_to_dictionary(B)}')

# Calculating pairings for different points
pair = pairing(A, B)
pair_AB = pairing(A, 2*B)
pair_BA = pairing(2*A, B)
pair_powered = pair**2

# We expect e(A,2B)=e(A,B)^2=e(2A,B)
print(f'\n\ne(A,B) = {pair}')
print(f'\n\ne(A,2B) = {pair_AB}')
print(f'\n\ne(2A,B) = {pair_BA}')
print(f'\n\ne(A,B)^2 = {pair_powered}')
