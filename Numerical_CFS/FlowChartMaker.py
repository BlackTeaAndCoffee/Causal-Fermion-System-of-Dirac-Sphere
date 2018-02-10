from pycallgraph import PyCallGraph
from pycallgraph import Config
from pycallgraph.output import GraphvizOutput
from pycallgraph import GlobbingFilter
from . import SymEngineFast
import symengine as si
import os
if( __name__ == "__main__"):
    config = Config(max_depth = 4)
    config.trace_filter = GlobbingFilter(include=[
            'C_F_S*','__main__','<module>','MainProg','ctypes*'
               ])
    graphviz = GraphvizOutput(output_file='PyCallGraphMinimizer.png')

    with PyCallGraph(output=graphviz, config=config):
        SymEngineFast.MainProg()
