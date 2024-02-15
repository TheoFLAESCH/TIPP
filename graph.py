####------------------------------------------------ graph.py -------------------------------------------####
# This programme is taking data from ROOT files, analysing and plot them using PYROOT. ------------------------------------###
# It works using functions defined in functions.py file. ###-------------------------------------------------------------------------------###
# It is required that the data files (.root) that we are going to analyse are in the same repertory than analysis.py and 
# functions.py. -----------------------------------------------------------------------### 
# One graph is displaying for 1 run (by default number 33) , 1 fitting method used (linear by default or sigmoid), 1 percent (by default 100%). These arguments are needed for running the script. -------------------------------------###

import argparse
from functions import *

parser = argparse.ArgumentParser(description='Plot graphs')
parser.add_argument('--numberRUN', default='33', help='Number of the run', required=False)
parser.add_argument('--percent', default='1', type=float, help='Percent of 12C ions (in percent). Examples : 0.2 for 20% of max(12C ions)', required=False)
parser.add_argument('--fittingMethod', default='linear', help='Fall-off adjustement method : linear or sigmoid ', required=False)
args = parser.parse_args()

##variables
numberRUN = args.numberRUN    #string
percent = args.percent    #string
fittingMethod = args.fittingMethod    #string


##list of files
data_file = [[],[]]   #array of strings
for i in os.listdir(os.getcwd()):
	if i.split("_")[0]=="ResultPB":
		data_file[0].append(i.split("_")[1].split(".")[0])
	if i.split("_")[0]=="TTreeVtx":
		data_file[1].append(i.split("_")[1])
	else:
		pass

###----------------------------------------------   MAIN  ---------------------------------------###
### FOR 1 RUN, 1 PERCENT OF 12C IONS, 1 FITTING METHOD (EITHER LINEAR OF SIGMOID) ##################


##get data from ROOT files
df_PB = GetFile(numberRUN, '/ResultPB_',data_file)
df_Vtx = GetFile(numberRUN, '/TTreeVtx_',data_file)

##fluctuation of number12Cions
f12C = df_Vtx.Histo1D(("f12C","f12C",500,0,16000000), "f12Cions")
z_nb_C, ordonn_nb_C  = GetData(f12C)[0], GetData(f12C)[2]
nb_12CMAX = z_nb_C[(np.array(ordonn_nb_C)).argmin()]  ##to get the max of the x axis of f12C histo
df_Vtx = df_Vtx.Filter( " f12Cions < %f*%f " %(nb_12CMAX,percent))
histo_vtx = df_Vtx.Histo1D(("histo_vtx","histo_vtx",500,-200,200),"fvtxZ")  ##number12C ions different from 100%
histo_vtx_cumule = histo_vtx.GetCumulative() # Cumulative histogram
prof1DPB = df_PB.Profile1D(("prof1D", "profile deposited energy", 500,0,200), "posz", "edep")

## Recovering data from histogrames in liste (easier to manipulate)

posz, err_posz, nombre_vtx, err_vtx = GetData(histo_vtx,80)
poszPB, err_poszPB, nombre_PB, err_PB = GetData(prof1DPB)
posz_c, err_posz_c, nombre_vtx_c, err_vtx_c = GetData(histo_vtx_cumule)

##PLOT from the same file. Creating TGraphErrors from these data and making them look good with TGraph_vtx_or_vtxc() function

g_vtx = ROOT.TGraphErrors(len(posz), posz, nombre_vtx, err_posz, err_vtx)
g_vtx_c = ROOT.TGraphErrors(len(posz_c), posz_c, nombre_vtx_c, err_posz_c, err_vtx_c)
g_PB = ROOT.TGraphErrors(len(poszPB), poszPB, nombre_PB, err_poszPB, err_PB)

TGraph_looking(g_vtx) # Grey curve
TGraph_looking(g_vtx_c) # Grey curve
TGraph_looking(g_PB) # Blue curve but cannot put it blue using the function (no idea why, we have to put it blue with the next two lines)


g_vtx.SetTitle("  ; Position (mm); Secondary vertex")
g_vtx_c.SetTitle("  ; Position (mm); Secondary vertex cumulated")
g_PB.SetLineColor(ROOT.kBlue)    #blue curve
g_PB.SetMarkerColor(ROOT.kBlue)

ROOT.gStyle.Reset("Modern")

print("______________Start of script______________\n")
print("Fitting of the fall-off of the vertices distribution\n")

print("--------------------------\n")	
print("{} fitting chosen\n".format(fittingMethod))
print("--------------------------\n")

if fittingMethod == 'linear':
	##linear fit part
	g_fit, g_y_over2 , z_inf, z_inf_err, a, b, delta_a, delta_b, chi_2_reduced = Fit_vtx(g_vtx, nombre_vtx, posz, 'L', 'graph')

if fittingMethod == 'sigmoid':	
	#sigmoid fit part
	g_fit, z_inf, z_inf_err, a, b, c, delta_a, delta_b, delta_c, chi_2_reduced = Fit_vtx(g_vtx, nombre_vtx, posz, 'S', 'graph')

if fittingMethod == 'cumulative':
	g_fit, z_inf, z_inf_err, a, b, c, delta_a, delta_b, delta_c, chi_2_reduced = Fit_vtx(histo_vtx_cumule, nombre_vtx, posz, 'C', 'graph')

print("--------------------------\n")
print("Chi_2_reduced = {} \n".format(chi_2_reduced))
print("Inflection point position = {} mm with error = {} \n".format(z_inf,z_inf_err))
print("Bragg peak position = {} mm \n".format(GetPB(prof1DPB)))

c = ROOT.TCanvas("c", "c", 1200, 800)
p1 = ROOT.TPad("p1", "", 0, 0, 1, 1)
p1.SetGrid()

p2 = ROOT.TPad("p2", "", 0, 0, 1, 1)
p2.SetFillStyle(4000) # will be transparent
p1.Draw()
p1.cd()

if fittingMethod == 'cumulative' :

	g_vtx_c.Draw("APL")
	g_vtx_c.GetHistogram().GetXaxis().SetTitleOffset(1.25)
	g_vtx_c.GetHistogram().GetYaxis().SetTitleOffset(1.25)
	ROOT.gPad.Update()

	tfont = g_vtx_c.GetHistogram().GetYaxis().GetTitleFont()
	tsize = g_vtx_c.GetHistogram().GetYaxis().GetTitleSize()
	lfont = g_vtx_c.GetHistogram().GetYaxis().GetLabelFont()
	lsize = g_vtx_c.GetHistogram().GetYaxis().GetLabelSize()

	leg = ROOT.TLegend(0.12, 0.5, 0.39, 0.7)
	leg.SetFillColor(ROOT.gPad.GetFillColor())
	leg.SetTextFont(lfont)
	leg.SetTextSize(lsize - 0.3 * lsize)
	leg.SetHeader("{} run ({} fit)".format(numberRUN,fittingMethod))
	leg.AddEntry(g_vtx_c, "Sec. vertex. cumul", "L")
	leg.AddEntry(g_PB, "Edep", "L")
	leg.Draw()

else :

	g_vtx.Draw("APL")
	g_vtx.GetHistogram().GetXaxis().SetTitleOffset(1.25)
	g_vtx.GetHistogram().GetYaxis().SetTitleOffset(1.25)
	ROOT.gPad.Update()

	tfont = g_vtx.GetHistogram().GetYaxis().GetTitleFont()
	tsize = g_vtx.GetHistogram().GetYaxis().GetTitleSize()
	lfont = g_vtx.GetHistogram().GetYaxis().GetLabelFont()
	lsize = g_vtx.GetHistogram().GetYaxis().GetLabelSize()

	leg = ROOT.TLegend(0.12, 0.5, 0.39, 0.7)
	leg.SetFillColor(ROOT.gPad.GetFillColor())
	leg.SetTextFont(lfont)
	leg.SetTextSize(lsize - 0.3 * lsize)
	leg.SetHeader("{} run ({} fit)".format(numberRUN,fittingMethod))
	leg.AddEntry(g_vtx, "Sec. vertex", "L")
	leg.AddEntry(g_PB, "Edep", "L")
	leg.Draw()

if fittingMethod == 'linear' : 
	leg.AddEntry(g_y_over2, "Half of Sec. \n vertex maximum", "L")
	g_y_over2.Draw()

g_fit.Draw()
ROOT.gPad.Update()
xmin = p1.GetUxmin()
xmax = p1.GetUxmax()
p1.cd()
dx = (xmax - xmin) / 0.8 #10 percent margins left and right
ymin = g_PB.GetHistogram().GetMinimum()
ymax = g_PB.GetHistogram().GetMaximum()
dy = (ymax - ymin) / 0.8  #10 percent margins top and bottom
p2.Range(xmin-0.1*dx, ymin-0.1*dy, xmax+0.1*dx, ymax+0.1*dy)
p2.Draw()
p2.cd()
g_PB.Draw("LP")
ROOT.gPad.Update()


axis = ROOT.TGaxis(xmax, ymin, xmax, ymax, ymin, ymax, 510, "+L")
axis.SetTitle("Deposit energy")   
axis.SetTitleFont(tfont)
axis.SetTitleSize(tsize)
axis.SetTitleColor(ROOT.kBlue)
axis.SetTitleOffset(1.25)
axis.SetLabelFont(lfont)
axis.SetLabelSize(lsize)
axis.SetLabelColor(ROOT.kBlue)
axis.SetLineColor(ROOT.kBlue)
axis.Draw()
ROOT.gPad.Update()
 
  
c.cd()
#c.Print("c{}.pdf".format(numberRUN))



print("--------------------- End of script ---------------------------------------------------------------------\n")

input()






