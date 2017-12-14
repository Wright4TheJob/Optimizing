# Optimizing
This repository hosts a variety of small optimization programs written in Python 3.

# Truss Optimization

The truss optimization program may be run by running the TrussMain.py file using Python 3.5 or newer. The GUI is built on PyQt5. This program minimizes the weight of a truss while not exceeding a specified stress in any member. Nodes may be selected as fixed or mobile for the optimization routine. Optimization of various cross section members coming soon.

# Four Bar Mechanism Synthesis

The four bar mechanism synthesis program will create a four bar mechanism (crank-rocker style currently) using a set of specified target points. The synthesis will match three points exactly and minimize error with the others. The program currently is performing symbolic differentiation of a large system of equations and is not viable for actual optimization. Future revisions will streamline this process significantly.
