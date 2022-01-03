from math import *

# parametry dla grs 80
a = 6378137
e2 = 0.00669437999013

b = a * (1 - e2)**0.5
f = 1 / 298.257222101

def M(phi):
    return (a*(1 - e2)) / ((1 - e2*sin(radians(phi))**2)**3)**0.5

def N(phi):
    return a / (1 - e2*sin(radians(phi))**2)**0.5

def pole(lam1, phi1, lam2, phi2):
    lam1 = radians(lam1)
    lam2 = radians(lam2)

    phi1 = radians(phi1)
    phi2 = radians(phi2)
    e = e2**0.5
    return  abs(
            (b**2 * (lam2 - lam1) / 2) * (((sin(phi2) /
            (1 - e2 * (sin(phi2)**2))) + (1 / (2 * e))
            * log(
                    (1 + e * sin(phi2)) / (1 - e * sin(phi2))))
                    - ((sin(phi1) / (1 - e2 * (sin(phi1)**2))) + (1 / (2 * e)) *
                log(
                    (1 + e * sin(phi1)) / (1 - e * sin(phi1)))))
            )

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

class Kivioj():
    def __init__(self, S, Az_p, phi_p, lam_p):
        ds = 1250
        n = int(S / ds)
        ds_o = S % ds
        for i in range(n):
            d_phi = (ds * cos(Az_p)) / M(phi_p)
            d_Az = (sin(Az_p) * tan(radians(phi_p)) * ds) / N(phi_p)

            Sr_phi = phi_p + degrees(d_phi) / 2
            Sr_az = Az_p + d_Az / 2
            d_phi_2 = (ds * cos(Sr_az)) / M(Sr_phi)
            d_lam = (ds * sin(Sr_az)) / (N(Sr_phi) * cos(radians(Sr_phi)))

            d_Az2 = (sin(Sr_az) * tan(radians(Sr_phi)) * ds) / N(Sr_phi)

            phi_p = phi_p + degrees(d_phi_2)
            lam_p = lam_p + degrees(d_lam)
            Az_p = Az_p + d_Az2

        d_phi = (ds_o * cos(Az_p)) / M(phi_p)
        d_Az = (sin(Az_p) * tan(radians(phi_p)) * ds_o) / N(phi_p)
        Sr_phi = phi_p + degrees(d_phi / 2)

        Sr_az = Az_p + d_Az / 2
        d_phi_2 = (ds_o * cos(Sr_az)) / M(Sr_phi)
        d_lam = (ds_o * sin(Sr_az)) / (N(Sr_phi) * cos(radians(Sr_phi)))

        d_Az2 = (sin(Sr_az) * tan(radians(Sr_phi)) * ds_o) / N(Sr_phi)
        phi_k = phi_p + degrees(d_phi_2)
        lam_k = lam_p + degrees(d_lam)

        Az = Az_p + d_Az2

        self.result = {
            "phi": phi_k,
            "lambda": lam_k,
            "Az": Az
        }


class Vincent():

    def __init__(self, lambda1, lambda2, fi_1, fi_2):
        U1 = atan((1 - f) * tan(radians(fi_1)))
        U2 = atan((1 - f) * tan(radians(fi_2)))

        L = radians(lambda2 - lambda1)

        difference = 100
        precision = radians(0.000001 / 3600)

        while difference >= precision:
            sins = ((cos(U2) * sin(L)) ** 2 + (
                    cos(U1) * sin(U2) - sin(U1) *
                    cos(U2) * cos(L)) ** 2) ** 0.5

            coss = sin(U1) * sin(U2) + cos(U1) * cos(U2) * cos(L)
            sig = atan(sins / coss)

            sina = (cos(U1) * cos(U2) * sin(L)) / sins
            cos2a = 1 - sina ** 2

            cos2sm = coss - ((2 * sin(U1) * sin(U2)) / cos2a)
            C = (f / 16) * cos2a * (4 + f * (4 - 3 * cos2a))

            prev_lambda = radians(lambda2 - lambda1) + (1 - C) * f * sina * (
                    sig + C * sins *
                    (cos2sm + C * coss * ((-1) +
                                          2 * cos2sm ** 2)))

            difference = abs(prev_lambda - L)
            L = prev_lambda

        u2 = ((a ** 2 - b ** 2) / b ** 2) * cos2a

        A = 1 + (u2 / 16384) * (4096 + u2 *
                                ((-768) + u2 * (320 - 175 * u2)))

        B = (u2 / 1024) * (256 + u2 *
                           ((-128) + u2 * (74 - 47 * u2)))

        d_sig = B * sins * (cos2sm + (1 / 4) * B * (
                coss * ((-1) + 2 * cos2sm ** 2) - (1 / 6) * B * cos2sm * ((-3) + 4 * sins ** 2) * (
                (-3) + 4 * cos2sm ** 2)))

        # odległość
        s = b * A * (sig - d_sig)
        self.s = s

        # azymuty
        y = cos(U2) * sin(prev_lambda)
        x = cos(U1) * sin(U2) - sin(U1) * cos(U2) * cos(prev_lambda)


        # poprawa azymutów
        if (y > 0 and x > 0):
            Az12 = atan(y / x)
        elif (y > 0 and x < 0):
            Az12 = atan(y / x) + pi
        elif (y < 0 and x < 0):
            Az12 = atan(y / x) + pi
        elif (y < 0 and x > 0):
            Az12 = atan(y / x) + 2 * pi

        self.Az12 = Az12
        y = cos(U1) * sin(prev_lambda)
        x = (-sin(U1)) * cos(U2) + cos(U1) * sin(U2) * cos(prev_lambda)
        # poprawa azymutó dla odwrotnego
        if (y > 0 and x > 0):
            Az21 = atan(y / x) + pi
        elif (y > 0 and x < 0):
            Az21 = atan(y / x) + 2 * pi
        elif (y < 0 and x < 0):
            Az21 = atan(y / x) + 2 * pi
        elif (y < 0 and x > 0):
            Az21 = atan(y / x) + 3 * pi

        self.Az21 = Az21





if __name__ == "__main__":
    a_fi = 50.25
    a_lam = 20.75
    b_fi = 50.25
    b_lam = 21.25
    c_fi = 50.
    c_lam = 20.75
    d_fi = 50.
    d_lam = 21.25



    a_d = Vincent(a_lam,d_lam, a_fi, d_fi)
    print("odleglosc AD i Azymut wprost: ",round(a_d.s,3), 'm', ',',  to_decimal(degrees(a_d.Az12)))
    middle_fi = (a_fi + d_fi) /2
    middle_lambda = (a_lam + d_lam) /2


    print("Punkt średniej szerokości: " +
          str((middle_fi)) + "°"
          + str((middle_lambda)) + "°"
          )
    # punkt środkowy

    kiv = Kivioj(a_d.s/2, a_d.Az12, a_fi, a_lam)
    print("Punkt środkowy AD: "
          + str((to_decimal(kiv.result['phi'])))
          + str(to_decimal(kiv.result['lambda']))

          )
    # pole powieszchni
    print("Pole powierzchni czworokata: "
          + str(round(pole(a_lam, a_fi, d_lam, d_fi),4)), 'm^2')


    srod_fi = kiv.result['phi']
    srod_lam = kiv.result['lambda']


    roznica = Vincent(middle_lambda,srod_lam, middle_fi, srod_fi)
    print("różnica odległości między punktem środkowym a średniej szerokości", round(roznica.s,3), 'm')
    print("Azymut wprost i odwrotny punktu średniej szerokości i punktu środkowego",
          to_decimal(degrees(roznica.Az12)) , "° ", to_decimal(degrees(roznica.Az21)), "°")


    #def __init__(self, lambda1, lambda2, fi_1, fi_2):
    # test = Vincent(20.75,20.75, 50.25, 50)
    # print(test.s)



