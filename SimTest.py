from pyx import *
import random
import numpy as np


def boltzmann(f_x, f_y, t, K):
    print('f_x,f_y', f_x, f_y)
    diff = abs(f_y - f_x)
    print('diff', diff)
    wert = np.exp(- diff/(t*K))
    #print 'f_x, f_y, diff, t, wert = ', f_x, f_y, diff, t, wert
    return wert

def nachbar(x_fitn, nachbaranzahl):
    nachbaranzahlhalb = int(nachbaranzahl/2)
    neg_maybe  = x_fitn - nachbaranzahlhalb
    much_maybe =  x_fitn  + nachbaranzahlhalb
    if neg_maybe < 0:
        ##print 1
        return nachbare1
    elif much_maybe > Stuetzstellen:
        ##print 2
        return nachbare2
    else:
        ##print 3
        #print x_int_list[x_fitn - nachbaranzahlhalb : x_fitn +nachbaranzahlhalb + 1]
        ##print 'x_fitn = ', x_fitn
        return x_int_list[x_fitn - nachbaranzahlhalb : x_fitn +nachbaranzahlhalb + 1]

def abschwingendeFunktion(tempus, rho, w, Amplitude):
    return np.exp(-rho*tempus**2)*(Amplitude*np.cos(w*tempus*0.5) + 1)

def f_gegn_x(x_values, y_values, X_Title, Y_Title, x_red, f_x, Quadratgroesse):

    x_min = min(x_values)
    x_max = max(x_values)
    y_min = min(y_values)
    y_max = max(y_values)

    c  = graph.graphxy(width=Quadratgroesse,height=Quadratgroesse,
                      x = graph.axis.linear(min = x_min, max = x_max, title= X_Title),
                      y = graph.axis.linear(min = y_min, max = y_max, title= Y_Title))
    c.plot(graph.data.values(x = x_values, y = y_values),[graph.style.line([color.rgb.blue])])
    faktor1 = Quadratgroesse/(x_max - x_min)
    faktor2 = Quadratgroesse/(y_max - y_min)

    c.fill(path.circle((x_red - x_min)*faktor1,(f_x - y_min)*faktor2, 0.1), [color.rgb.red])
    return c

def temperatur(Iter, hilfsarray_fuer_temp, rho, w, Amplitude):
    temp =[]
    for k in range(Iter):
        temp.append(abschwingendeFunktion(hilfsarray_fuer_temp[k],rho,w,
            Amplitude))
    return temp

def Bilder(hilfsarray_fuer_temp, temp, berg_tal, x_fitn, fitn_wert_x, m, s):
    leinwand = canvas.canvas()

    leinwand2 = f_gegn_x(x, berg_tal, 'r=radius', 'landschaft', x[x_fitn]
            , fitn_wert_x,15)
    leinwand3 = f_gegn_x(hilfsarray_fuer_temp, temp, 'zeit=t', 'temperatur',
    hilfsarray_fuer_temp[m], temp[m], 3)
    leinwand.insert(leinwand2)
    leinwand.insert(leinwand3,[trafo.translate(10,12)])
    leinwand.writeGSfile('./BilderZuSImuli/Simuliannili%04d'%s, device = 'png16m')

def Bild_berechnetes_Minimum(x_min):

    leinwand = canvas.canvas()
    leinwand3 = f_gegn_x(x, berg_tal, 'radius = r', 'Landschaft',
    x[x_min[0]], x_min[1], 10)
    leinwand.insert(leinwand3)
    leinwand.writePDFfile('Loesung')


def Minimum_finden():
    x_fitn = random.randint(0, Stuetzstellen - 1) #Startwert... random
    fitn_wert_x = berg_tal[x_fitn]
    x_min = [x_fitn, fitn_wert_x]
    for m,t in enumerate(temp):
        for j in range(10):
            randomi = random.random()
            fitn_wert_y = berg_tal[nachbar(x_fitn, nachbaranzahl)[random.randint(0,nachbaranzahl - 1)]]
            boltzi = boltzmann(fitn_wert_x, fitn_wert_y, t)
            if fitn_wert_x > fitn_wert_y:
                x_fitn = berg_tal.index(fitn_wert_y)
                if fitn_wert_y < x_min[1]:
                    x_min[0]=x_fitn
                    x_min[1]=fitn_wert_y
                fitn_wert_x=fitn_wert_y

            elif (fitn_wert_x < fitn_wert_y) and (randomi < boltzi) :
                x_fitn = berg_tal.index(fitn_wert_y)
                if fitn_wert_y < x_min[1]:
                    x_min[0]=x_fitn
                    x_min[1]=fitn_wert_y

                fitn_wert_x = fitn_wert_y
    print(x_y_min[0], x_y_min[1])


if __name__ == "__main__":
    '''
    Dieses Programm dient der Einarbeitung in das Simulated Annealing Verfahren.
    Alles andere ist selbst erklaerend. ;)
    '''

    '''
    Stuetzstellen fuer die Funktion deren Minimum gefunden werden soll.
    '''

    Stuetzstellen=300

    x= list(np.linspace(0,2*np.pi,Stuetzstellen))  #x-Achse fuers Diagramm
    x_int_list=[i for i in range(Stuetzstellen)]   #x_achse fuer die Rechnung, hier koennen
                                                   #direkt mit integer Zahlen rechnen

    Iter = 5                                      #Die Anzahl der temperaturiterationen
    hilfsarray_fuer_temp = np.linspace(0.01,25,Iter)

    Amplitude = 0.5                                #Amplitude der Fluktuation
                                                   #der Temperatur
    w = np.pi                                    #Die Frequenz mit der die temp variieren soll
    rho = 0.001                                      #Der Halbwertswert fuer die abschwingfunktion

    w_b1 = 3.14                                       #Die erste Frequenz der
    w_b2 = 7.17                                       #die zweite Frequ. d. Lands.

    '''
    Das Vorgehen ist folgendermassen. Es wird ueber die verschiedenen Temperaturen iteriert.
    Anfangen sollte man mit hoeheren Temperaturen. Was hohe Temperaturen sind, sollten in diesem
    Zusammenhang nochmal untersucht werden. Je nachdem wie hoch die Huegel oder Temperaturen sind
    und welche Nachbarn alles beruecksichtigt werden, muessten vielleicht  hoehere Temperaturen beruecksichtigt
    werden. Bei jeder Temperatur werden dann nochmal mehrere Iterationen durchgefuehrt. Hier gibt es
    ja nur diskrete Zeitpunkt und nichts kontinuierliches, daher sollte man bei jeder Temperatur
    eine gewisse Anzahl an Iterationen durchfuehren. Alles andere ist selbsterklaerend.
    '''
    berg_tal = berg_Und_Tal(Stuetzstellen, w_b1, w_b2)
    temp = temperatur(Iter, hilfsarray_fuer_temp, rho, w, Amplitude)
    nachbaranzahl = Stuetzstellen

    nachbare1 = [i for i in range(nachbaranzahl)] #nachbar fuer den Fall, dass der momentane
                                                  #x-wert zu nahe an Null ist... verstehst scho

    hilfspoint =Stuetzstellen - nachbaranzahl     #
    nachbare2 = [hilfspoint  + i for i in range(nachbaranzahl)]


    Minimum_finden(Bild = False)
    #Minimum_finden_und_Bilder()
