import pkg_resources, os


def get_data(path):
	ressource_name = os.path.join("data", path)
	return pkg_resources.resource_filename("Numerical_CFS", ressource_name)

def get_config(path):
	ressource_name = os.path.join("config", path)
	return pkg_resources.resource_filename("Numerical_CFS", ressource_name)
