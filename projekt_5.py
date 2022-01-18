from math import *
import numpy as np

# grs 80
e2 = 0.00669437999013
a = 6378137
b = a * (1 - e2) ** 0.5
f = 1 / 298.257222101

# krasowski
akr = 6378245
bkr = 6356863
e2kr = (akr ** 2 - bkr ** 2) / (akr ** 2)

def set_m(phi):
    return (a * (1 - e2)) / ((1 - e2 * sin(radians(phi)) ** 2) ** 3) ** 0.5


def set_n(phi):
    return a / (1 - e2 * sin(radians(phi)) ** 2) ** 0.5




def to_x_y_z(fi, lam, h, A, e2):  # podac w radianach
    N = A / (1 - (e2) * sin(fi) ** 2) ** 0.5

    x = (N + h) * cos(fi) * cos(lam)
    y = (N + h) * cos(fi) * sin(lam)
    z = ((N * (1 - e2) + h) * sin(fi))
    return round(x, 3), round(y, 3), round(z, 3)





def Hirvonen(x, y, z):
    # krasowski
    akr = 6378245
    bkr = 6356863
    e2kr = (akr ** 2 - bkr ** 2) / (akr ** 2)

    r = ((x ** 2) + (y ** 2)) ** 0.5

    # first iteration
    first_fi = atan(z / r * (1 - e2kr) ** -1)  # radiany

    N = akr / (1 - e2kr * sin(first_fi) ** 2) ** 0.5
    h = r / cos(first_fi) - N


    second_fi = atan(z / r * (1 - e2kr * N / (N + h)) ** -1)

    epsilon = 2.42406840554768e-10


    while abs(second_fi - first_fi) >= epsilon:
        N = akr / (1 - e2kr * sin(first_fi) ** 2) ** 0.5
        h = r / cos(first_fi) - N

        first_fi = second_fi
        second_fi = atan(z / r * (1 - e2kr * N / (N + h)) ** -1)


    lam = atan(y / x)
    N = akr / (1 - e2kr * sin(second_fi) ** 2) ** 0.5

    h = r / cos(second_fi) - N
    #print("kontrola h", r, cos(second_fi), N), h
    #check_x, check_y, check_z = to_x_y_z(second_fi, lam, h, akr, e2kr)

    #print("kontrola danych dla hiveriona, orginał: ",x, y, z)
    #print("ponownie przeliczone: ", round(check_x,2), round(check_y,2), round(check_z,2))
    # print(check_x, check_y, check_z)

    return to_decimal(degrees(second_fi)), to_decimal(degrees(lam)), round(h, 3)  # degrees, degrees


def BursyWolf(xp, yp, zp):
    kappa = 0.8407728e-6
    alfa = -1.786877784465417e-06
    beta = -2.5612706773016787e-07
    gam = 4.089597325631379e-06

    x0 = -33.4297
    y0 = 146.5746
    z0 = 76.2865

    first = np.array([xp, yp, zp])
    second = np.array(
        [[kappa, gam, -beta],
         [-gam, kappa, alfa],
         [beta, -alfa, kappa]
         ])
    third = np.array([x0, y0, z0])

    ret = first + second @ first + third

    for i in range(0, len(ret)):
        ret[i] = round(ret[i], 3)
    return ret

def to_decimal(decimal_degree):
    stopnie = int(decimal_degree)
    minuty = int((decimal_degree-stopnie) * 60)
    sekundy = (decimal_degree - stopnie - minuty/60)*3600

    sekundy = round(sekundy,5)
    sekundy = str(sekundy)
    stopnie = str(stopnie)
    minuty = str(minuty)

    # korekta zer
    if len(stopnie) == 1:
        stopnie = "0" + stopnie

    if len(minuty) == 1:
        minuty = "0" + minuty

    if len(sekundy) == 1:
        sekundy = "0" + sekundy

    return stopnie + '° ' + minuty + "' " + sekundy + "'' "

if __name__ == '__main__':
    A = {'fi': 50.25, 'lam': 20.75}
    B = {'fi': 50.00, 'lam': 20.75}
    C = {'fi': 50.25, 'lam': 21.25}
    D = {'fi': 50.00, 'lam': 21.25}

    E = {'fi': 50.125, 'lam': 21.00}
    F = {'fi': 50.12527044824481, 'lam': 21.000651088258433}
    points = [A, B, C, D, E, F]

    # print(Hirvonen(10,20,30))
    # print(BursyWolf(10,20,30))
    for P in points:

        print('punkt:', P)
        print("dla grs 80")
        x, y, z = to_x_y_z(radians(P['fi']), radians(P['lam']), 100, a, e2)  # default h = 100 m
        print(x, y, z)
        print("dla elipsoidy krasowskiego")
        xk, yk, zk = BursyWolf(x, y, z)
        print(xk, yk, zk)
        print('po ponownym przeliczeniu na fi, lambda, h')
        print(Hirvonen(xk, yk, zk))

