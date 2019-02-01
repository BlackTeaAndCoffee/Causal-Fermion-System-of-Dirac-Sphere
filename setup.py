from setuptools import  setup, find_packages

setup(
    name = "Num_CFS_Cluster",
    version = "0.2",
    description = "Basically for a set of parameters i use Symengine to produce a function, which is the integrand of the Action-Integral. The integration is done numerically. I am searching for the set of parameters, which minimizes the Integral.",
    author = "Mustafa Akman",
    author_email = "", # FIXME: add email
    license = "GPL v3",
    install_requires = ["numpy", "sympy", "scipy"],
    package_data = {"Numerical_CFS": ["config/*", "data/*"]
        },
    packages = find_packages())

