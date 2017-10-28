import os

_root = os.path.abspath(os.path.dirname(__file__))
def get_data(path):
	return os.path.join(_root, "data", path)

def get_config(path):
	return os.path.join(_root, "config", path)
