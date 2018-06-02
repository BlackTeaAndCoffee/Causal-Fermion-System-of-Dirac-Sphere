from setuptools import  setup, find_packages

setup(
    name = "Numerical_CFS_Docker",
    version = "0.1.9.8.1",
    description = "Basically for a set of parameters i use Symengine to produce a function, which is the integrand of the Action-Integral. The integration is done numerically. I am searching for the set of parameters, which minimizes the Integral.",
    author = "Mustafa Akman",
    author_email = "mustisummer@gmail.com", 
    license = "GPL v3",
    install_requires = ["numpy", "symengine", "scipy"],
    package_data = {"Numerical_CFS": ["config/*", "data/*"]
        },
    packages = find_packages())

