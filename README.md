# Causal-Fermion-System-of-Dirac-Sphere
Basically for a set of parameters i use Symengine to produce a function, which is the integrand of the Action-Integral. The integration is done numerically. I am searching for the set of parameters, which minimizes the Integral.

[Causal-Fermion-Systems](https://en.wikipedia.org/wiki/Causal_fermion_system) are yet another approach to describe 
fundamental physics. 
This code is being developed for my Phd-Thesis. It is here for all who are interested in my work and 
maybe want to see how i did the calculations, use it for their own work or give me any useful hints.
This work is based on the Master Thesis of Nikki Kilbertus, of which you can download a copy on the wiki-page of this project.

The Main Programm is SymEngineFast.py. I Called it like this, because i use [SymEngine](https://github.com/symengine), which is a computer algebra package, which is faster than [Sympy](http://www.sympy.org/en/index.html) in contructing the integrand. The reason for me in constructing the integrand algebraicly is
- Reproduction of Niki's Results
- if it's fast enough this way, it's awesome, because the numerical error only appear in the integration. 
  So i only need to read the error approximation given by the integration routine. 
- If i construct the whole integrand numericaly, currently i have no 
  good idea on how to approximate the error. 
     - So that would be a good issue to work on!

For the integration i at first used the quadpack routine of scipy. The symbolic integrand needed to get lambdifiyd for that. 
But all of this takes to much time. So i went over to generate C-Code with the print.ccode function of sympy and then used ctypes.scipy to integrate, which is all in all faster. 
- The generating of the ccode takes for larger sets a lot of time, maybe there is a way to that faster.
- I never went really into discovering, which integration routine is the most sutiable, because i was very 
fixated on finding out, if this is at all feasable. 

The second Part of the Code would be the Minimizing-Problem. For this i currently use the Simulated-Annealing-Algorithm. But this has some down sides, because for that i need to calculate the action for a given set of parameters a lot of times. The action is the integral of the integrand over a fixed space. 
    
   - The integration takes depending on the variables too much time. I never wrote benchmarks of that, so that would be another issue. 
   - Another Problem is that for larger sets of paramters, the amount of needed integrations also increases. 
 
There are a lot of different minimizing methods out there, so another one would be bayesian optimization. It's completely different from Simulated-Annealing and maybe a good alternative. 
# Required Installations  
 The Main Programms are in Numerical_CFS. 
 Install following:
 - Python3 (Numpy, Scipy)
 - Sympy
 - SymEngine 
 - I use [Pyx](http://pyx.sourceforge.net/) for plotting, but you are free to choose any other package for plotting. 
 Pyx is at first not so easy to use, but the moment you got used to it, you can do
                                        a lot of things with it. 
                                        
 - You need to be able to run C programms. 
  
 - To use the FlowChartMaker.py you need to install [PyCallgraph](http://pycallgraph.slowchop.com/en/master/#)
 

# Pip
For Linux-Users, type into Terminal:
`pip install Numerical_CFS`
# or 
Download the package "tar.gz" file from [Numerical_CFS](https://pypi.python.org/pypi/Numerical-CFS) unpack it and then 
run `python setup.py install`


