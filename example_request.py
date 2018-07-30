# coding: utf-8
import requests
data = {"Vary_Parameters": (True, True, False), 
"Set_Initialstate_randomly": (False, False, False),
"Impuls": (0, 15, "[iF for iF in range(1)]"),
"Frequenz": (0, 20, "[iF for iF in range(1)]"),
"Constraints": (1, 0.005, "[1]"),
"System_sizes": (1, 1, 1),
"Test": (True, False)}
r = requests.post('http://127.0.0.1:8080/calculate', json=data)
