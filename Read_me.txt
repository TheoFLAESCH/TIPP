
####################### QUICK EXPLENATION OF THE SUBJECT #################################

Hadrontherapy is possible because charged particle lose most of their energy at a given distance, this feature is called the Bragg Peak.
Along their path in the body, they also interact and produce secondary particule which are detectable and it is possible to reconstruct their vertex of interaction.
The vertices distribution is rapidly increasing when the beam enters the material, then it reaches a plateau (small or long, depending on the beam energy) and decreases rapidly, generally 
close to the Bragg peak. Our objectif is to find a relationship between the Bragg peak and the inflection point of the fall-off region of the vertex distribution.

####################### FILE "graph.py" #################################

This is probably the first code to run to understand what we have done and the methods

This code has for objective to plot the data that are given to us which are : The energy deposition as a function of the position (gives the Bragg peak) and the vertex distribution 
(from which we want to recover the inflection point of the fall-off region).

This code plot on a same graph these two data, with the left Y axis the number of vertex, the right Y axis the energy deposition and in X axis the depth in the body along the beam axis.
In addition it will plot the fit of the fall-off region depending on the chosen method, and it is also possible to select a certain run (energy deposition + vertex for a beam energy).

If you execute the file with "python graph.py" by default it will use the run 33 data and use the linear fitting method.

If you want to change the run, add the argument --numberRUN followed by the number of the run. Example : python graph.py --numberRUN 34
There are 6 runs possible : 33, 34, 37, 38, 39, 42

If you want to change the function fitting the fall-off region add the argument --fittingMethod followed by the function. Example : python graph.py --fittingMethod linear
There are 2 possible function to fit the vertex : linear, sigmoid.

There is a third method, but it uses the cumulative of the vertex distribution. As said in the report, the fit doesn't converge, however you can still plot for any run this method
to see what it looks like using the argument --fittingMethod cumulative

If you want to change the number of vertices, add the argument --percent. For example : python graph.py --percent 0.2 (to take 20% of the vertices)

You can naturally use the 3 arguments at once. For example : python graph.py --numberRUN 37 --fittingMethod sigmoid --percent 0.8


/!\ Taking a percentage of vertex too low may result in not enough data for the fit to converge making the code longer to run.
/!\ The code may run for a bit of time with the sigmoid and cumulative function because the fit doesn't converge for all run as mentionned in the report.


This code also return in the terminal the reduced chi square of the fit (the chi square divided by the number of degres of freedom), the position of the Bragg peak, the estimated inflection
point and the uncertainties.

We remind you that the goal is to run over all the set of data corresponding to different beam energy {energy loss, vertex} and to recover {inflection point, Bragg peak} to find a 
relationship between the inflection point and the Bragg peak, which a next file will do.

####################### FILE "function.py" #################################

This is the file in which all the functions used for graph.py are stored in.

