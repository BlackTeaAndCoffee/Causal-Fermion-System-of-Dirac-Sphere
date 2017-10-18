from pycallgraph import PyCallGraph
from pycallgraph import Config
from pycallgraph.output import GraphvizOutput
from pycallgraph import GlobbingFilter
from SymEngineFast import MainProg
import os

config = Config()
config.trace_filter = GlobbingFilter(include=[
        'SymEngineFast*'
            ])
graphviz = GraphvizOutput(output_file='PyCallGraph.png')

with PyCallGraph(output=graphviz, config=config):
    MainProg()
