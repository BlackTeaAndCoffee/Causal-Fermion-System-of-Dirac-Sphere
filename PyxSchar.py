from pyx import *

def initialXY(xmin, xmax, ymin, ymax, x_title, y_title,  Height, Width, keypos):
    return graph.graphxy(width=Width,height=Height,
        x = graph.axis.linear(min = xmin, max = xmax,
                          title = x_title ),
        y = graph.axis.linear(min = ymin,max = ymax,
                          title= y_title),
                      key = graph.key.key(pos = keypos, dist =0.1))


def set_values(x_list, y_list, Kurve_Name):
    return graph.data.values(x = x_list, y =  y_list, title = Kurve_Name)


def plot_diag(ListData, Title, TitlePos, PDFName, canvas):
    canvas.plot(ListData,[graph.style.line([color.gradient.Rainbow])])
    canvas.text(*TitlePos, Title)
    canvas.writePDFfile(PDFName)

