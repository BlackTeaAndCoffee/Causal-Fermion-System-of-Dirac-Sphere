from SymEngineErweiterung import get_Wirkung
import Diag_Schar
import Rho_data
import symengine as si
import PyxSchar
import numpy as np


def func(t, r, N, Intgrenze, T, Rho_Liste,w_Liste,kappa,k1):
    K_Liste[ind]= k1
    print(K_Liste, WertListe, 'yaaaaaaay1')
    Wert = get_Wirkung(t, r, N, Intgrenze, T, K_Liste, Rho_Liste,w_Liste,kappa)
    WertListe.append(Wert[0])
    KK.append(k1)
    print(WertListe, 'yaay2')
    return WertListe, KK

if __name__ == "__main__":
    si.var('r t')
    T = 1 #Lebensdauer des Universums, wird fuer die Schwartzfunktion benoetigt
    N = 2

    kappa = 0.01
    kappa_Anzahl = 1

    x_Anf = 0
    x_End = np.pi
    Intgrenze = [x_Anf, x_End]

    ind = 0 #Welcher Impuls soll variert werden, 0 für K1 und 1 für K2


    w_Liste = [0.,1.]
    K_Liste = [0.,0.01]
    AnzahlGewichte = 2
    Kurve_Names = []
    PDFTitle = '$GewichteScharVarK%2.2fs$'%ind

    c = PyxSchar.initialXY(0,20,0,1.5, r'$K_2$','Wirkung',10,10 , 'tr')

    dd=[]

    for k in range(AnzahlGewichte):
        rho1 = 0.1*k
        Rho_Liste = [rho1, (1- rho1)/3]
        dk = 0.1
        k1 = 0.
        WertListe = []
        KK =[]
        while 1:
            WertListe, KK = func(t, r, N, Intgrenze, T, Rho_Liste,
                w_Liste,kappa, k1)
            l1 = len(WertListe) - 1
            if  l1>=2 :
                a = abs((WertListe[l1-2] - WertListe[l1-1])/(KK[l1-2] -
                    KK[l1-1]))
                b = abs((WertListe[l1] - WertListe[l1-1])/(KK[l1] - KK[l1-1]))
                if  b != 0 and (a/b < 0.5 or a/b >2):
                    WertListe.pop(l1)
                    KK.pop(l1)
                    k1 -= dk
                    dk = 0.02
                    k1 += dk
                    WertListe, KK = func(t, r, N, Intgrenze, T, Rho_Liste,
                w_Liste,kappa, k1)
                elif a != 0 and (b/a < 0.5 or b/a > 2):
                    WertListe.pop(l1)
                    KK.pop(l1)
                    k1 -= dk
                    dk = 0.02
                    k1 += dk
                    WertListe, KK = func(t, r, N, Intgrenze, T, Rho_Liste,
                w_Liste,kappa,  k1)
                else:
                    dk = dk*2
                    dk = min(dk, 0.5)
            if WertListe[l1] > 1.5 or k1 > 20:
                break

            k1 += dk
        LL = np.array(WertListe)
        LL = LL/LL[0]
        print(LL, KK)
        dd.append(PyxSchar.set_values(KK, LL, r"$\rho_1$=%2.2f"%rho1))


    if ind ==0:
        PyxSchar.plot_diag(dd, '$kappa=%2.5f, K_%2.2f=%2.2f$'%(kappa, 2, 0),
            (5,10), PDFTitle, c)
    else:
        PyxSchar.plot_diag(dd, '$kappa=%2.5f, K_%2.2f=%2.2f$'%(kappa, 1, 0),
            (5,10), PDFTitle, c)

