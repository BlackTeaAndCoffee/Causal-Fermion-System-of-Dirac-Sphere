from pyx import *


def Plot(NameOfFile, NameOfPDF):
    g = graph.graphxyz(size=4, projector=graph.graphxyz.parallel(50, 25),
            x=graph.axis.lin(title=r"$K_1$"),
            y=graph.axis.lin(title=r"$K_2$"),
            z =graph.axis.lin(title = r"Wirkung"))
    g.plot(graph.data.file(NameOfFile, x=1, y=2, z=3, color=3),
           [graph.style.surface(gradient=color.gradient.RedGreen,
                                gridcolor=color.rgb.black,
                                backcolor=color.rgb.black)])
    g.writePDFfile(NameOfPDF)

if __name__ == "__main__":
    N  = 2
    rho1 = 0.9
    Rho_Liste = [rho1, (1-rho1)/3]
    Plot("foo.txt", "WirkunN%dRho1_%f_Rho2_%f_VarK1_K2" %(N, Rho_Liste[0], Rho_Liste[1]))

