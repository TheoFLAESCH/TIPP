###----------------------------------- functions.py FILE-------------------------------------------###
###-------------------------- File containing all functions needed --------------------------------###

import ROOT 
import os
import sys
from array import array
import copy
import numpy as np


###-------------------------------------------------------------------------------###
#return df from file depending on char= ResultPB_ or TTreeVtx_
def GetFile(nbRUN,char,data_file):
	if char == "/ResultPB_" :
		if [x==nbRUN for x in data_file[0]]==[False]*len(data_file[0]):
			print("No file")	
			sys.exit(0)
		else:
			data_file = os.getcwd() + char + nbRUN + ".root"
			df = ROOT.RDataFrame("Element", data_file)
	if char == "/TTreeVtx_" :
		if [x==nbRUN for x in data_file[1]]==[False]*len(data_file[1]):
			print("No file")
			sys.exit(0)
		else:
			data_file = os.getcwd() + char + nbRUN + "_1_Norm.root"
			df = ROOT.RDataFrame("tree", data_file) 
	else :
		pass
	
	return df

###-------------------------------------------------------------------------------###
#get data from histo (ordon=ordonnee,z=bins, shift for having the same origin for the two graphs)
def GetData(histo, shift=0):
	z, err_z = array("d"), array("d")
	ordon, err_ordon = array("d"), array("d")
	for i in range(1, histo.GetNbinsX()):
		z.append(histo.GetBinCenter(i) + shift)
		err_z.append(1)
		ordon.append(histo.GetBinContent(i))
		err_ordon.append(histo.GetBinError(i))
	err_z = array("d", [(max(z)-min(z))/histo.GetNbinsX()])*len(z)
	return z, err_z,ordon, err_ordon

###-------------------------------------------------------------------------------###
#Get maximum value of vertex distribution by averaging (used in linear and sigmoid fittings)
def GetMaxVtx(nb_vtx, threshold = 0.95):
	y_max = max(nb_vtx)
	y_tsd = threshold*y_max
	above_y_tsd = array("d")
	for i in nb_vtx:
		if i > y_tsd:
			above_y_tsd.append(i)
		else:
			pass
	y_max_av = 0   #y_max average
	for i in above_y_tsd:
		y_max_av = y_max_av + i
	y_max_av = y_max_av/len(above_y_tsd)
	
	return y_max_av

###-------------------------------------------------------------------------------###
##return index in x of the closest value in liste 
#list = array("d"), and K = float
def returnZclosest(liste,K):
	list_copy = copy.deepcopy(liste)
	R = np.abs((np.array(liste)-GetMaxVtx(liste))).argmin()
	for i in range(0,R):    ##to focus only on the right side of the curve
		list_copy.pop(0)
	tab = np.abs(np.array(list_copy)-K)
	idx = (tab).argmin()
	return R + idx 

###-------------------------------------------------------------------------------###
##Get peak bragg, return value of the bragg peak position
def GetPB(histo):
	binmax = histo.GetMaximumBin() #more populated bin
	z = histo.GetXaxis().GetBinCenter(binmax) #the corresponding abscisse
	return z 


###-------------------------------------------------------------------------------###
##Avoid repetition in main. 
def TGraph_looking(g):
	if g == 'g_vtx' or 'g_vtx_c': #grey curve
		#if g == 'g_vtx' :
		#	g.SetTitle("  ; Position (mm); Secondary vertex")
		#if g == 'g_vtx_c' :
		#	g.SetTitle("  ; Position (mm); Secondary vertex cumulated")
		
		g.SetMarkerStyle(8)
		g.SetMarkerSize(0.75)
		g.SetMarkerColor(16)
		g.SetLineColor(16)   #grey curve
		g.SetFillColor(3)

	if g == 'g_PB': #blue curve
		g.SetMarkerStyle(8)
		g.SetMarkerSize(0.75)
		g.SetFillColor(3)

###----------------------------- METHOD 1 : Linear fit --------------------------###
###----------------------------- METHOD 2 : Sigmoid fit -------------------------###
###----------------------------- METHOD 1 : Linear fit --------------------------###
###----------------------------- METHOD 2 : Sigmoid fit -------------------------###
###----------------------------- METHOD 3 : Cumulative fit ----------------------###
def Fit_vtx(g_vtx, nombre_vtx, posz, method, return_graph):
	y_max = GetMaxVtx(nombre_vtx)  #keep the same way to get y_max
	y1 = y_max * 3 / 10   #25% of y_max
	y2 = y_max * 7 / 10 #75% of y_max
	z1 = posz[returnZclosest(nombre_vtx,y1)]
	z2 = posz[returnZclosest(nombre_vtx,y2)]

	y1c = y_max * 25 / 100   #30% of y_max
	y2c = y_max* 75 / 10 #70% of y_max
	z1c = posz[returnZclosest(nombre_vtx,y1c)]
	z2c = posz[returnZclosest(nombre_vtx,y2c)]

	if method == "L" :

		f1 = ROOT.TF1("f1", '[0]* x + [1]', z1, z2)
		f1.SetLineWidth(8)  
		g_vtx.Fit('f1', 'RN')    #"N" fit not drawn

		#Get Parameters
		a = f1.GetParameter(0)
		b = f1.GetParameter(1)
		chi_2_reduced = f1.GetChisquare() / f1.GetNDF()
		delta_a = f1.GetParError(0)
		delta_b = f1.GetParError(1)
		Y = (y_max/2) - b 
		z_inf = (1 / a) * Y  #Inflection point  (found analytically)
		z_inf_err = np.sqrt(Y**2*delta_a**2/a**4  +  delta_b**2/a**2)
		
		if return_graph == 'graph' :

			## to extend the linear fit
			x, y, y_over2 = array("d"), array("d") , array("d") 
			for i in range(int(z2)-20,int(z1)+20):
				x.append(i)
				y.append(a*i+b)
				y_over2.append(y_max/2)     #y_max getting from GetMaxVtx()
	       
			g_fit = ROOT.TGraph(len(x), x, y)
			g_fit.SetMarkerStyle(8)
			g_fit.SetMarkerSize(0.75)
			g_fit.SetMarkerColor(ROOT.kRed)
			g_fit.SetLineColor(ROOT.kRed)    #red curve
			g_fit.SetFillColor(3)
			g_fit.SetLineWidth(3)

			g_y_over2 = ROOT.TGraph(len(x), x, y_over2)
			g_y_over2.SetMarkerStyle(8)
			g_y_over2.SetMarkerSize(0.75)
			g_y_over2.SetMarkerColor(109)
			g_y_over2.SetLineColor(8)    #green curve
			g_y_over2.SetFillColor(3)
			g_y_over2.SetLineWidth(3)

			return g_fit,g_y_over2, z_inf, z_inf_err, a, b, delta_a, delta_b, chi_2_reduced
		else :
			return z_inf, z_inf_err, a, b, delta_a, delta_b, chi_2_reduced

	if method == "S" : 

		f1 = ROOT.TF1("f1", '[0] / ([1] + exp(- [2] *x))', z1, z2)

		f1.SetParameter(0,y_max)
		f1.SetParameter(1,0)
		f1.SetParameter(2,-0.01)

		f1.SetLineWidth(8)  
		g_vtx.Fit('f1', 'RN')    #"N" fit not drawn

		#Get Parameters
		a = f1.GetParameter(0)
		b = f1.GetParameter(1)
		c = f1.GetParameter(2)
		chi_2_reduced = f1.GetChisquare() / f1.GetNDF()
		delta_a = f1.GetParError(0)
		delta_b = f1.GetParError(1)
		delta_c = f1.GetParError(2)
		Y = (y_max/2) - b 
		z_inf = - np.log(b) / c   #Inflection point (found analytically)
		z_inf_err = np.sqrt((delta_b**2)/(b**2*c**2)  +  (z_inf**2*delta_c**2)/(c**2))
		
		if return_graph == 'graph' :

			## to extend the sigmoid fit
			x, y = array("d"), array("d") 
			for i in range(int(z2) ,int(z1)):
				x.append(i)
				y.append(a / (b + np.exp(- c * i)))

		       	g_fit = ROOT.TGraph(len(x), x, y)
			g_fit.SetMarkerStyle(8)
			g_fit.SetMarkerSize(0.75)
			g_fit.SetMarkerColor(ROOT.kRed)
			g_fit.SetLineColor(ROOT.kRed)    #red curve
			g_fit.SetFillColor(3)
			g_fit.SetLineWidth(3)

			return g_fit, z_inf, z_inf_err, a, b, c, delta_a, delta_b, delta_c, chi_2_reduced  
		else :
			return z_inf, z_inf_err, a, b, c, delta_a, delta_b, delta_c, chi_2_reduced

	if method == "C" :
		
		f1 = ROOT.TF1("f1", '([0] * log([1] * exp([2] * x) + 1)) / ([1] * [2]) + [3]', z1c, z2c)

		f1.SetParameter(0, y_max)
		f1.SetParameter(1, 0.1)
		f1.SetParameter(2, - 0.01)
		f1.SetParameter(3, g_vtx.GetMaximumBin())

		g_vtx.Fit('f1', 'RN') #Do not draw

		a = f1.GetParameter(0)
		b = f1.GetParameter(1)
		c = f1.GetParameter(2)
		d = f1.GetParameter(3)
		chi_2_reduced = f1.GetChisquare() / f1.GetNDF()
		delta_a = f1.GetParError(0)
		delta_b = f1.GetParError(1)
		delta_c = f1.GetParError(2)
		delta_d = f1.GetParError(3)
		z_inf = - np.log(b) / c   #Inflection point (found analytically)
		z_inf_err = np.sqrt((delta_b**2)/(b**2*c**2)  +  (z_inf**2*delta_c**2)/(c**2))
		
		if return_graph == 'graph' :

			## to extend the sigmoid fit
			x, y = array("d"), array("d") 
			for i in range(int(z2) ,int(z1)):
				x.append(i)
				y.append((a * np.log(b * np.exp(c * i) + 1)) / (b * c) + d)

		       
			g_fit = ROOT.TGraph(len(x), x, y)
			g_fit.SetMarkerStyle(8)
			g_fit.SetMarkerSize(0.75)
			g_fit.SetMarkerColor(ROOT.kRed)
			g_fit.SetLineColor(ROOT.kRed)    #red curve
			g_fit.SetFillColor(3)
			g_fit.SetLineWidth(3)

			return g_fit, z_inf, z_inf_err, a, b, c, delta_a, delta_b, delta_c, chi_2_reduced  
		else :
			return z_inf, z_inf_err, a, b, c, delta_a, delta_b, delta_c, chi_2_reduced








