#!/usr/bin/env python
import ROOT
import sys
import math
import numpy
import ctypes
import argparse

#check presence of a given key in root file 
def check_isInFile(name, rootfile):
	isInFile_key=True
	key = rootfile.FindKey(name)
	if (key==None) :
		print '\n',name,'does not exists in ',  rootfile, '\n'
		isInFile_key=False
	return isInFile_key

def get_BinsX(graphname, whichBins):
	npoints = graphname.GetN()
	#print 'npoints :', npoints
	bin_s =[]

	if (whichBins == 'xlow'):
			for iPoint in range(0,npoints): 
				bin_s.append(graphname.GetErrorXlow(iPoint))	
	if (whichBins == 'xhigh'):
			for iPoint in range(0,npoints): 
				bin_s.append(graphname.GetErrorXhigh(iPoint))	
	if (whichBins == 'xnominal'):
			for iPoint in range(0,npoints): 
				bin_s.append(graphname.GetX()[iPoint])
	return bin_s


def get_BinsY(graphname, whichBins):
	npoints = graphname.GetN()
	#print 'npoints :', npoints
	bin_s =[]

	if (whichBins == 'ylow'):
			for iPoint in range(0,npoints): 
				bin_s.append(graphname.GetErrorYlow(iPoint))	
	if (whichBins == 'yhigh'):
			for iPoint in range(0,npoints): 
				bin_s.append(graphname.GetErrorYhigh(iPoint))	
	if (whichBins == 'ynominal'):
			for iPoint in range(0,npoints): 
				bin_s.append(graphname.GetY()[iPoint])
	return bin_s


def get_BinsXlowEdges(graph):
	Nbins = graph.GetN()
	if (Nbins<1):
		sys.exit('graph'+graph+' has zero bins!')
	err_low = get_BinsX(graph, 'xlow')
	nominal = get_BinsX(graph, 'xnominal')
	lowEdges= [center - error for center, error in zip(nominal,err_low)]
	return lowEdges

def get_BinsXupEdges(graph):
	Nbins = graph.GetN()
	if (Nbins<1):
		sys.exit('graph'+graph+' has zero bins!')
	err_up = get_BinsX(graph, 'xhigh')
	nominal = get_BinsX(graph, 'xnominal')
	upEdges= [center + error for center, error in zip(nominal,err_up)]
	return upEdges
		
def have_sameBinning (graph1, graph2):
	N_1 = graph1.GetN()
	N_2 = graph2.GetN()
	if (N_1 != N_2):
		print 'Can not proceed! Graphs do not have same number of points.'
		print 'Points in', graph1, ': ', N_1, 'Points in', graph2, ': ', N_2
		return False
	
	if ( len(set(get_BinsXlowEdges(graph1))-set(get_BinsXlowEdges(graph2)))==0 and len(set(get_BinsXupEdges(graph1))-set(get_BinsXupEdges(graph2))) ==0 ):
		return True
	else:
		return False


def make_ratio(data, mc):
	ratio_values = []
	for iBin in range(0,data.GetN()):
		if (mc.GetY()[iBin]==0.):
			ratio_values.append(1000)
			print 'bin ', iBin , 'in ', mc, 'has content = 0'
		else:  
			ratio_values.append((data.GetY()[iBin])/(mc.GetY()[iBin]))	
	return ratio_values 

def average_error(graph,point): #per point
	e_high = graph.GetErrorYhigh(point)
	e_low = graph.GetErrorYlow(point)
	return ((e_high+e_low)/2)

def calculate_ratioError (num, den): 
	error_s = []
	for iPoint in range(0,num.GetN()):
		sigmaN = average_error(num, iPoint)
		sigmaD = average_error(den, iPoint)
		error2 = 0.
		error2 = ((sigmaN/num.GetY()[iPoint])**2)
		error2 = error2 + ((num.GetY()[iPoint]*sigmaD/(den.GetY()[iPoint]**2))**2)
		error_s.append(math.sqrt(error2))
	return error_s
		
		
#create canvas with the ratio graph
def plotRatio(title, xtitle, ytitle, plotName, format_s, binsX, errXlow, errXhigh, binsY, errY):
	can = ROOT.TCanvas(title, title, 800, 600)
	N = len(binsX)
	# converting python lists to c arrays of double* 
	# needed for the tgraph constructor
	BINSx = (ctypes.c_double * len(binsX))(*binsX)
	BINSy = (ctypes.c_double * len(binsY))(*binsY)
	ERRxLOW = (ctypes.c_double * len(errXlow))(*errXlow)
	ERRxHIGH =  (ctypes.c_double * len(errXhigh))(*errXhigh)
	ERRy = (ctypes.c_double * len(errY))(*errY)
	ratio_graph = ROOT.TGraphAsymmErrors(N, BINSx, BINSy, ERRxLOW, ERRxHIGH, ERRy, ERRy)
	#style settings
	ratio_graph.SetMarkerColor(4)
	ratio_graph.SetMarkerStyle(21)
	ratio_graph.SetMarkerSize(1.2)
	ratio_graph.SetTitle(title)
	ratio_graph.GetXaxis().SetTitle(xtitle)
	ratio_graph.GetXaxis().SetTitleOffset(1.2)
	ratio_graph.GetYaxis().SetTitle(ytitle)
	ratio_graph.GetYaxis().SetTitleOffset(1.3)
	can.cd()
	ratio_graph.Draw('ALP')
	#save canvas 
	for aFormat in format_s:
		can.SaveAs(plotName+'.'+aFormat)
	

def printText(what, variable, what_ratio, textFileName, binsX, low_edges, up_edges, data_values, data_errHigh, data_errLow, mc_values, mc_errLow, mc_errHigh, ratio_values, ratio_errors):


	f = open(textFileName+'.txt','w') 
	f.write('\n'+what+' as a function of '+variable+'\n\n')
	f.write('-------------------------------------------------------------------------------------------------------------------------\n')
	f.write(variable+'|\tScaleFactor\t\t|\t\tData\t\t\t|\t\tMC\t\t\t|\n' )
	f.write('-------------------------------------------------------------------------------------------------------------------------\n')
	for i in range(0,len(binsX)):
		f.write( str(low_edges[i])+' - '+str(up_edges[i])+' |\t') # pt bin  |
		f.write("{0:.4f}".format(ratio_values[i]) )# scale factor + - err |
		f.write(' + - '+"{0:.4f}".format(ratio_errors[i]) )
		f.write('\t|\t') 
		f.write("{0:.4f}".format(data_values[i]) ) # data         + err - err |
		f.write(' + '+"{0:.4f}".format(data_errHigh[i]) )
		f.write(' - '+"{0:.4f}".format(data_errLow[i])  )
		f.write('\t|\t') 
		f.write("{0:.4f}".format(mc_values[i]) ) # data         + err - err |
		f.write(' + '+"{0:.4f}".format(mc_errHigh[i]) )
		f.write(' - '+"{0:.4f}".format(mc_errLow[i])  )
		f.write('\t|\t') 
		f.write('\n')
	f.write('-------------------------------------------------------------------------------------------------------------------------\n')
	f.close()




#######################################################################


#################################################
'''
arguments to be put in a parser (now hard coded)
	- inputRootFile 
	- Graph Name
	- what (e.g. Electron ID Efficiency)
	- variable (e.g. electron pT)
	- makePlot (True, False)
	- formats ( e.g. ['png', 'root'] )
	- makePlainText (True,False)
'''
#################################################

parser = argparse.ArgumentParser()

parser.add_argument("inputRootFile", help="full name of input root file (e.g. filename.root)")
parser.add_argument("TGraphName", help="name of the TGgraphs in the input root file. Assuming it contains name_Data and name_MC")

'''
parser.add_argument("--what", action="store", dest="what", default="efficiency", help="what is contained in the input file (e.g.: Electron ID Efficiency)")

parser.add_argument("--ratio", action="store", dest="what_ratio", default="ratio",help="what is the ratio going to be (e.g. : Scale Factor Data/MC) ")
parser.add_argument("--var", action="store", dest="variable", default="x", help="what is on the x axis (e.g. : Electron pT)")

parser.add_argument('--plot', action="store", dest="makePlot", type=bool, default=True, help=" save the plot (True, False)")

parser.add_argument("--plotName", action="store", nargs='+', default="plot", dest="plotName", help="name of the plot") 

parser.add_argument("--formats", nargs='+', action="store", type=str, default="png", dest="format_s",help = "formats to save the plot (e.g.: png root)") 

parser.add_argument("--txt", type=bool, action="store", default="False", dest="makePlainText", help="save a plain text file with graph contents (True, False)")
parser.add_argument("--tname", action="store", default="table", dest="textFileName",help="name of the output text file")#change to optional
'''
args = parser.parse_args()


print 'args.what_ratio', args.what_ratio
print 'args.variable', args.variable
 

##################################################
# names 

what = 'Electron ID Efficiency'
what_ratio = 'Scale Factor Data/MC'
variable = 'Electron p_T'
makePlot = False
format_s = ['png', 'root']
plotName = 'plot'
makePlainText = True
textFileName = 'values'

#################################################

inRootFile = ROOT.TFile(args.inputRootFile, 'read')

suffix = ["_Data", "_MC"]
#formats for plotting

format_s = ["png","root"]

# check existence of graphs
for aSuffix in suffix:
	keyname = args.TGraphName+aSuffix
	foundGraphs = check_isInFile(keyname, inRootFile)
	if not foundGraphs:
		sys.exit('could not find '+args.TGraphName+' in root file' )

#read graphs 
graph_data = inRootFile.Get(args.TGraphName+suffix[0])
graph_mc = inRootFile.Get(args.TGraphName+suffix[1])

#check that they have same binning 
if (have_sameBinning(graph_data,graph_mc)==False):
	sys.exit('Can not proceed. Graphs do not have same binning. Please check')

#make the ratio 
ratio_values = make_ratio(graph_data, graph_mc) #list containing bin by bin ratio

#calculate error on the ratio 
ratio_errors = calculate_ratioError(graph_data, graph_mc)

#save useful quantities
binsX = get_BinsX(graph_data, 'xnominal')
errXlow = get_BinsX(graph_data, 'xlow')
errXhigh = get_BinsX(graph_data, 'xhigh')
up_edges = get_BinsXupEdges(graph_data)
low_edges =  get_BinsXlowEdges(graph_data)
data_values  = get_BinsY(graph_data, 'ynominal')
data_errHigh = get_BinsY(graph_data, 'yhigh')
data_errLow  = get_BinsY(graph_data, 'ylow') 
mc_values    = get_BinsY(graph_mc,   'ynominal')
mc_errHigh   = get_BinsY(graph_mc,   'yhigh')
mc_errLow    = get_BinsY(graph_mc,   'ylow')

#draw scale factors in a canvas 
if args.makePlot == True:
	if (len(format_s)>0) :
		plotRatio(what, variable, what_ratio, plotName, format_s, binsX, errXlow, errXhigh, ratio_values, ratio_errors)

#store values in a plain text file 
if args.makePlainText == True:
	printText(what, variable, what_ratio, textFileName, binsX, low_edges, up_edges, data_values, data_errHigh, data_errLow, mc_values, mc_errLow, mc_errHigh, ratio_values, ratio_errors)



