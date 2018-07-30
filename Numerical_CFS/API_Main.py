import cherrypy
from io import StringIO
from .Lib_Action_Minimum import MainProg

class API(object):
    def __init__(self):
        pass
    @cherrypy.expose
    @cherrypy.tools.json_in()
    def calculate(self):
        def configfunktion(name):
            return cherrypy.request.json[name]

        file = StringIO()
        
        MainProg(configfunktion, file)
        return file


if( __name__ == "__main__"):
    cherrypy.quickstart(API())

