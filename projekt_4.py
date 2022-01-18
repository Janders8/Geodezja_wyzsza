from math import *
from shapely.geometry import Polygon

#grs 80
e2 = 0.00669437999013
e2p = 0.00673949674226495
a = 6378137
b = a * (1 - e2)**0.5
#f = 1 / 298.257222101
L0_1992 = L0=19*pi/180 # radians


def L0_for_2000(nr):
    if nr == 8:
        L0 = radians(24)
    elif nr == 7:
        L0 = radians(21)
    elif nr == 6:
        L0 = radians(18)
    elif nr == 5:
        L0 = radians(15)

    return L0

def set_m(phi): # podac jako radiany
    return (a*(1 - e2)) / ((1 - e2*sin(phi)**2)**3)**0.5

def set_n(phi): # podac jako radiany
    return a / (1 - e2*sin(phi)**2)**0.5


def xyGK(F, La, uklad):

    result = [0,0]
    if uklad == '1992':
        L0=19*pi/180
    elif uklad == '2000':
        if La >= 22.5:
            nr = 8
            L0 = radians(24)
        elif La <= 22.5 and La > 19.5:
            nr = 7
            L0 = radians(21)
        elif La <= 19.5 and La > 16.5:
            nr = 6
            L0 = radians(18)
        elif La <16.5:
            nr = 5
            L0 = radians(15)
        result.append(nr)

    #print(L0)
    F= radians(F)
    La = radians(La)
    N = set_n(F)

    t=tan(F)
    mi2=(e2p)*((cos(F))**2)

    A0=1-(e2/4)-((3*(e2)**2)/64)-((5*(e2)**3)/256)
    A2=(3/8)*((e2+((e2**2)/4)+((15*(e2)**3))/128))
    A4=(15/256)*(e2**2+((3*((e2)**3))/4))
    A6=(35*((e2)**3))/3072

    sigma=a*(A0*F-A2*sin(2*F)+A4*sin(4*F)-A6*sin(6*F))
    l=La-L0

    Xgk=sigma+((l**2)/2)*N*sin(F)*cos(F)*(1+((l**2)/12)*
        ((cos(F))**2)*(5-(t**2)+9*(mi2)+4*((mi2)**2))+
        ((l**4)/360)*((cos(F))**4)*(61-58*(t**2)+(t**4)+
        270*(mi2)-330*(mi2)*(t**2)))

    Ygk=l*N*cos(F)*(1+((l**2)/6)*((cos(F))**2)*
            (1-(t**2)+(mi2))+((l**4)/120)*((cos(F))**4)*
            (5-18*(t**2)+(t**4)+14*(mi2)-58*(mi2)*(t**2)))

    result[0] = Xgk
    result[1] = Ygk

    return result

def to_1992(F, La):
    m0_1992=0.9993
    Xgk, Ygk = xyGK(F, La, '1992')

    X=m0_1992*Xgk-5300000
    Y=m0_1992*Ygk+500000


    return (X,Y)

def from_1992_to_gk(x, y):
    m0_1992=0.9993
    Xgk = (x +5300000)/ m0_1992
    Ygk = (y - 500000) / m0_1992

    return (Xgk, Ygk)


def from_2000_to_gk(x, y, nr):
    m0_2000=0.999923

    Xgk = x / m0_2000

    Ygk = (y-(nr*1000000) - 500000) / m0_2000


    return (Xgk, Ygk)



def to_2000(F, La):
    Xgk, Ygk, nr = xyGK(F, La, '2000')
    m0_2000=0.999923

    x=m0_2000*Xgk
    y=m0_2000*Ygk+1000000*nr+500000
    return (x,y,nr)



def gk_to_geo(x,y, L0 ):

    A0=1-(e2/4)-((3*(e2)**2)/64)-((5*(e2)**3)/256)
    A2=(3/8)*((e2+((e2**2)/4)+((15*(e2)**3))/128))
    A4=(15/256)*(e2**2+((3*((e2)**3))/4))
    A6=(35*((e2)**3))/3072

    fi_0 = x / (a * A0)

    eps = 100
    while eps > ((0.00001 / 3600) * (pi / 180)):
        fi = (fi_0 + (x - (a * ((A0 * fi_0) - (A2 * sin(2 * fi_0)) + (A4 * sin(4 * fi_0)) - (A6 * sin(6 * fi_0))))) / (a * A0))
        eps = abs(fi - fi_0)
        fi_0 = fi


    # t to tangens fi
    t = tan(fi)
    t2 = t**2

    mi2 = (e2p) * ((cos(fi)) ** 2)
    N = set_n(fi)
    M = set_m(fi)

    fi_final = fi - (((y ** 2) * t) / (2 * M * N))\
                * (1 - ((y ** 2) / (12 * (N ** 2)))
                * (5 + 3 * t2 + mi2 - 9 * mi2
                * t2 - 4 * (mi2 ** 2)) +
                ((y ** 4) / (360 * (N ** 4))) *
                ( 61 + 90 * t2 + 45 * (t2 ** 2)))

    lam_final = L0 + (y / (N * cos(fi))) * (1 - ((y**2) / (6*(N**2))) *
                (1 + 2*t2 + mi2) + ((y**4) / (120*(N**4)))
                * (5 + 28*t2 + 24*(t2**2) + 6*mi2 + 8*mi2*t2))

    return  degrees(fi_final), degrees(lam_final)

def skala_dlugosci_1992(ygk, fi): # to z tabelki zaokraglone # - > m układu , kappa w kilometrach

    fi = radians(fi)
    M = set_m(fi)
    N = set_n(fi)
    q = (M*N)**0.5


    mgk = 1 + ((ygk**2) / (2*q**2)) + ((ygk**2)/(24*q**4))


    m92 = 0.9993 * mgk

    #kappa
    Z = ( m92-1) *1000

    m_kwadrat = 0.9993**2 * mgk**2

    Z_ha = (m_kwadrat -1) * 10000

    return(round(m92,6), round(Z,3), round(m_kwadrat,6), round(Z_ha,3))

def skala_dlugosci_2000(ygk, fi):

    fi = radians(fi)
    M = set_m(fi)
    N = set_n(fi)
    q = (M*N)**0.5


    mgk = 1 + ((ygk**2) / (2*q**2)) + ((ygk**2)/(24*q**4))

    m2000 = 0.999923 * mgk

    Z = ( m2000-1) *1000


    m_kwadrat = 0.999923**2 * mgk**2

    Z_ha = (m_kwadrat -1) * 10000

    return(round(m2000,6), round(Z,3), round(m_kwadrat,6), round(Z_ha,3))

def skala_dlugosci_gk(ygk, fi):

    fi = radians(fi)
    M = set_m(fi)
    N = set_n(fi)
    q = (M*N)**0.5


    mgk = 1 + ((ygk**2) / (2*q**2)) + ((ygk**2)/(24*q**4))



    Z = (mgk -1) *1000


    m_kwadrat = mgk**2

    Z_ha = (m_kwadrat -1) * 10000

    return(round(mgk,6), round(Z,3), round(m_kwadrat,6), round(Z_ha,3))

if __name__ == '__main__':
    a_fi = 50.25
    a_lam = 20.75

    c_fi = 50.25
    c_lam = 21.25

    b_fi = 50.
    b_lam = 20.75

    d_fi = 50.
    d_lam = 21.25

    e_fi = 50.125
    e_lam = 21.0

    f_fi = 50.12527044824481
    f_lam = 21.000651088258433



    fis = [a_fi, b_fi, c_fi, d_fi, e_fi, f_fi ]
    lams = [a_lam, b_lam, c_lam, d_lam, e_lam, f_lam]

    x_1992 = []
    y_1992 = []

    x_2000 = []
    y_2000 = []
    y_2000_numbers = []

    Xgk = []
    Ygk = []



    # przeliczanie pkt

    for i in range(len(fis)):
        x, y = xyGK(fis[i],lams[i], '1992')
        Xgk.append(x)
        Ygk.append(y)

        f, l = to_1992(fis[i],lams[i])
        x_1992.append(f)
        y_1992.append(l)

        f, l, nr = to_2000(fis[i],lams[i])
        x_2000.append(f)
        y_2000.append(l)
        y_2000_numbers.append(nr)


    print("pkt")
    print("punkty w gausie-krugerze: ", )
    for i in range(len(fis)):
        print(round(Xgk[i],3), round(Ygk[i], 3))
    print("punkty w układzie 1992: ", )
    for i in range(len(fis)):
        print(round(x_1992[i],3), round(y_1992[i],3))

    print("punkty w układzie 2000: ", )
    for i in range(len(fis)):
        print(round(x_2000[i],3), round(y_2000[i],3))

    #pola
    print("pola")
    cords_1992 = [(x_1992[0], y_1992[0]), (x_1992[1], y_1992[1]), (x_1992[3], y_1992[3]), (x_1992[2], y_1992[2]),
                  (x_1992[0], y_1992[0])]
    cords_2000 = [(x_2000[0], y_2000[0]), (x_2000[1], y_2000[1]), (x_2000[3], y_2000[3]), (x_2000[2], y_2000[2]),
                  (x_2000[0], y_2000[0])]
    cords_gk = [xyGK(a_fi, a_lam, "1992"), xyGK(b_fi, b_lam, "1992"), xyGK(d_fi, d_lam, "1992"),
                xyGK(c_fi, c_lam, "1992"), xyGK(a_fi, a_lam, "1992")]


    #print("wsp poligonu xy: ",cords_1992)
    poly_1992 = Polygon(cords_1992)
    poly_2000 = Polygon(cords_2000)
    poly_gk = Polygon(cords_gk)
    print("pole na elipsoidzie gausa-grugera", poly_gk.area/1000000)
    print("pole w układzie 1992", poly_1992.area/1000000)
    print("pole w układzie 2000", poly_2000.area/1000000)

    # droga powrotna
    new_1992_fi =[]
    new_1992_lam = []
    new_2000_fi = []
    new_2000_lam = []

    #dla 1992
    new_x1992_gk = []
    new_y1992_gk = []

    new_x2000_gk = []
    new_y2000_gk = []

    for i in range(len(fis)):
        #1992
        xgk, ygk = from_1992_to_gk(round(x_1992[i],3), round(y_1992[i],3))
        new_x1992_gk.append(xgk)
        new_y1992_gk.append(ygk)
        fi, lam = gk_to_geo(xgk,ygk, L0_1992)
        new_1992_fi.append(fi)
        new_1992_lam.append(lam)
        #2000
        xgk, ygk = from_2000_to_gk(round(x_2000[i],3), round(y_2000[i],3), y_2000_numbers[i])
        new_x2000_gk.append(xgk)
        new_y2000_gk.append(ygk)
        fi, lam = gk_to_geo(xgk,ygk, L0_for_2000(y_2000_numbers[i]))
        new_2000_fi.append(fi)
        new_2000_lam.append(lam)

        # gk


        #print(new_1992_fi[i],new_1992_lam[i])


        print("skale zniekształceń")
        print("dla GK")
        for i in range(len(new_1992_fi)):
            print(skala_dlugosci_gk(Ygk[i], gk_to_geo(Xgk[i], Ygk[i], L0_1992)[0]))
            #
            #
        print("dla 1992")
        for i in range(len(new_1992_fi)):
            print(skala_dlugosci_1992(new_y1992_gk[i], new_1992_fi[i]))

        print("dla 2000")
        for i in range(len(new_1992_fi)):
            print(skala_dlugosci_2000(new_y2000_gk[i], new_2000_fi[i]))


