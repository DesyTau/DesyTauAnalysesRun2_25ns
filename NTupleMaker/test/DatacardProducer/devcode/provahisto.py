import ROOT
from ROOT import TH1D, TH2D

def unroll(histo2D):			
	if histo2D is not None:
		htitle = histo2D.GetTitle()
		hname = histo2D.GetName()
		n_bins_x = (histo2D.GetNbinsX())*(histo2D.GetNbinsY())
		histo1D = TH1D(hname, htitle, n_bins_x, 0, n_bins_x)	
		for i in xrange(1,n_bins_x+1):
			histo1D.SetBinContent(i, histo2D.GetBinContent(i))
			histo1D.SetBinError(i, histo2D.GetBinError(i))
		return histo1D
		#histo2D = histo1D #unrolled result
	else:
		pass


h = TH2D("h","h",3,0.,3.,2,0.,2.)
print "nbins_x : " , h.GetNbinsX() , "nbins_y: ", h.GetNbinsY()
h.SetBinContent(0,0)
#for i in xrange(1, (h.GetNbinsX()*h.GetNbinsY())+1):
#	print "i: ", i
#	h.SetBinContent(i,-i)

'''
h.SetBinContent(0.5,0.5,-1)
h.SetBinContent(1.5,0.5,-2)
h.SetBinContent(1.5,0.5,-3)
h.SetBinContent(1.5,0.5,-3)
h.SetBinContent(1.5,0.5,-3)
'''

for y in range(0,h.GetNbinsY()+2):
	for x in range(0,h.GetNbinsX()+2):
		print " x : ", x , " y: ", y , " - ", h.GetBin(x,y), " - ", h.GetBinContent(h.GetBin(x,y))

h1 = unroll(h)

h1.SaveAs("h1.root")
h.SaveAs("h.root")
