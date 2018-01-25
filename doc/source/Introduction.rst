Introduction
************

Currently i am working on the numerical Side of `Causal-Fermion-Systems <http://en.wikipedia.Causal_fermion_system>`_. To be precise currently all of this packages are for a specific class of Causal-Fermion-Systems. Namely the CFS of a Dirac-Sphere in 4 Dimensions. In Future i will hopefully add different Systems as well and then alter the code to do so accordingly. 

One big aspect of Causal-Ferimion-Systems is the Variation. But codewise here i basically for a set of parameters use Symengine to produce a function, which is the integrand of the Action-Integral. The integration is done numerically. I am searching for the set of parameters, which minimizes the Integral.

Conceptionally speaking this programm can be split into get_Integrand, get_Action, get_Minimum.
In get_Integrand as the name already spills it out, i'm constructing the
Integrand, which will be later on integrated in get_Action to get the Action. In get_Minimum
an algorithm finds ideally the minimum.

get_Integrand
=============
    Here my idea was to use sympy to get the analytic form of the
    integrand. So this way ideally only in the numeric integration
    i would get round off errors. The integration is done numerically
    because analytically it takes to much time. I have to do this integration
    a lot of times, for various parameters. So for smaller integrands it
    could be feasable to do the integration analytically but the integrand will
    get bigger and bigger for N getting  bigger. We are talking about several
    pages. So in comparison doing the integration numerically seems more
    plausible. I am not really able to say, how long the analytic integration
    would take for an expression, so i have no exlicit remarks on that.

    Another option would have been, to do everything with numpy arrays.
    I have done this once, if i remember correctly. This wasn't fast
    at all. The problem is, that i still have to construct the integrand
    out of huge amount of smaller terms. I think python as the overhead (i'm
    not sure if thats the correct name for it) causes it to be that slow.
    This code is in FullyNumeric.py.
    
    *You can find this function in SymEngineFast.py.*


get_Action
==========
    In get_Action i integrate over the Integrand i get via Get_Integrand.
    Depending on which method i decided to to use, different methods come
    to use.
    There is integrating with nquad from scipy, or nquad ctypes, or integrating
    with c.
    
    *You can find this function in SymEngineFast.py.*


get_Minimum
===========
    This Programm contains the code for minimizing the Action. With the
    help of some Minimum finder, i want to find ideally the global
    Minimum of the Action. The Action as encoded in SymEngineFast.py has
    the parameters K_List, Rho_List and w_List.  By parameters i mean,
    the elements of these Lists.

    The Paramater space is in general infinite dimensional. Our idea is
    to take N=1 and the corresponding sub-parameter-space and find the
    minimum there. Then we expand to the subspace corresponding to N=2,
    which contains the subspace for N=1 and find the Minimum there.
    Ideally we wish to find some convergency for higher getting N.

    The Minimizer we choose is the **simulated annealing procedure**.
    The reason for that is, that the Action must not be smooth.
	
    * A Problem with that, which already seems to get an issue is that
      each evaluation of the action takes to much time. On Average
      higher N means longer evaluation time. For higher N we also need
      to probe the Paramater space more often. So the simulated
      annealing method might not be good enough.

    * There is another method called bayesian optimization. I haven't
      tried this one out yet. But i will definitely do so soon.

    * In the simulated annealing algorithm we basically begin with some
      initial state. A state is just some combinations of the various
      paramters in K_List, Rho_List and w_List.  The initial state
      should either be given by the user or the programm should pick
      some state.
		 
           + The user would give a state, if he or she want's to give
             the algorithm a starting point. This we will do by each
             increase of the Shell-Number N.

    * A state always gets an Energyvalue (float number) and the state
      with the lowest Energy is the desired one. Now we can imagine,
      that the energy landscape has a lot of mountains and valleys and
      sometime you need to pass a huge mountain to get to the deepest
      valley. For this kind of problem the simulated annealing
      algorithm is very good. Usualy after a change of variables, the
      combination of parameters with the lowest energy gets picked as
      the state to be in, but if the temperature is high enough
      sometimes the algorithm would also jump to state in which the
      energy is higher.  For more details, please look "Simulated
      annealing" up.  
      
    *You can find this function in Action_Minimum.py*	
	
FYI
===
     For more information on these Functions, look in to the master
     thesis of Nikki Kilbertus, which is included on the
     github site. If you are up to it, you can read from page 1, but from a
     numerical perspective Chapter 2.3 Numerical Recipe (page 24ff) is important.


