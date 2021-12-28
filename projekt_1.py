import math as m
import numpy as np
import pandas as pd
import geopandas as gpd
import matplotlib.pyplot as plt

# współrzędne lotniska w seulu
s_fi, s_lambda, s_h = 37.453628, 126.423247, 10

e2 = 0.00669438002290
# grs80 a parameter
A = 6378137

F1 = np.deg2rad(s_fi)
L1 = np.deg2rad(s_lambda)

# macierz obrotu
R = np.array(
    [[-m.sin(F1) * m.cos(L1), -m.sin(L1), m.cos(F1) * m.cos(L1)],
     [-m.sin(F1) * m.sin(L1), m.cos(L1), m.cos(F1) * m.sin(L1)],
     [m.cos(F1), 0, m.sin(F1)]
     ])
R = R.transpose()

def to_x_y_z(fi, lam, h, A, e2):
    fi = np.deg2rad(fi)
    lam = np.deg2rad(lam)
    N = A / (1 - (e2) * m.sin(fi) ** 2) ** 0.5

    x = (N + h) * m.cos(fi) * m.cos(lam)
    y = (N + h) * m.cos(fi) * m.sin(lam)
    z = ((N * (1 - e2) + h) * m.sin(fi))
    return x, y, z

#           wsp lotniska wsp lotu
def geo2neu(F1, L1, H1, F2, L2, H2, R):
    pkt_1 = to_x_y_z(F1, L1, H1, A, e2)
    pkt_2 = to_x_y_z(F2, L2, H2, A, e2)


    x = np.array([[pkt_2[0] - pkt_1[0]],
                  [pkt_2[1] - pkt_1[1]],
                  [pkt_2[2] - pkt_1[2]]])
    neu = R @ x
    return neu


def to_azimuth(down, up, func):
    azi = func(up/down)
    if down > 0 and up > 0:
        azi= np.rad2deg(azi)
    elif down < 0 and up > 0:
        azi=np.rad2deg(azi + m.pi)
    elif down < 0 and up < 0:
        azi=np.rad2deg(azi + m.pi)
    else:
        azi= np.rad2deg(azi + 2*m.pi)
    if azi > 360:
        azi -=360
    elif azi < 0:
        azi += 360
    return azi






def main():
    # wczytanie danych
    lot = np.loadtxt("lot.txt")
    dane_lotu = pd.DataFrame(lot)


    #wyświetlenie trasy lotu we współrzędnych geodezyjnych fi lambda h
    airport_gdf = gpd.GeoDataFrame(dane_lotu, geometry = gpd.points_from_xy(dane_lotu[1], dane_lotu[0]))
    airport_gdf.plot(markersize = 1.5, figsize = (10, 10))
    plt.show()


    # konwersja na neu
    n = []
    e = []
    u = []

    for i in lot:
        temp = geo2neu(s_fi, s_lambda, s_h, i[0], i[1], i[2], R )
        n.append(temp[0][0])
        e.append(temp[1][0])
        u.append(temp[2][0])

        # wypisanie kolejnych współrzędnych neu
        #print(temp[0][0],temp[1][0],temp[2][0])

    # wyświetlenie neu
    fig = plt.figure()
    ax = plt.axes(projection="3d")
    ax.scatter3D(n, e, u, c=u, cmap='cividis')




    #wykres współrzędnej u
    #plt.plot(range(len(u)), u)

    # skośna, A i Z
    for i in range(1, len(n)):

        root = (n[i]**2 + e[i]**2 + u[i]**2)**0.5

        print(to_azimuth(n[i],e[i],np.arctan),
              root,
              to_azimuth(root,u[i] , np.arccos))


    print("^ Azymut A, odległość skośna oraz Kąt Z ^")






    for a, b, c in zip(n, e, u):
        if c < 0:
            ax.scatter(a, b, c, "red", s= 200 )
            ax.text(a, b, c, '%s' % ("       samolot znika za horyzontem"), size=10, zorder=1, color='k')
            print("samolot zniknął za horyzontem we wsp neu: ", a, b, c)
            break

    ax.text(n[-1], e[-1], u[-1], '%s' % ("     warszawa"), size=10, zorder=1, color='k')
    ax.text(n[0], e[0], u[0], '%s' % ("     seul "), size=10, zorder=1, color='k')
    plt.show()

if __name__ == "__main__":
    main()
