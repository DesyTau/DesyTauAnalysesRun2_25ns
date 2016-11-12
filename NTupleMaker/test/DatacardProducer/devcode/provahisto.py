import ROOT
from ROOT import TH1D, TH2D

def unroll(histo2D):			
	if histo2D is not None:
		htitle = histo2D.GetTitle()
		hname = histo2D.GetName()
		n_bins_x = (histo2D.GetNbinsX())*(histo2D.GetNbinsY())
		histo1D = TH1D(hname, htitle, n_bins_x, 0, n_bins_x)

		k = 1
		for j in xrange(1, histo2D.GetNbinsY()+1):
			for i in xrange(1, histo2D.GetNbinsX()+1):
				histo1D.SetBinContent(k, histo2D.GetBinContent(histo2D.GetBin(i,j)))
				histo1D.SetBinError(k, histo2D.GetBinError(histo2D.GetBin(i,j)))
				print "i: ", i, " j: ",j,"| Set content of bin ",k, " : ", histo2D.GetBinContent(histo2D.GetBin(i,j))
				k+=1	
		return histo1D
	else:
		pass


h = TH2D("h","h",3,0.,3.,2,0.,2.)
print "nbins_x : " , h.GetNbinsX() , "nbins_y: ", h.GetNbinsY()
#h.SetBinContent(0,0)
print "--content of 2D histogram --"

for j in xrange(1, h.GetNbinsY()+1):
	for i in xrange(1,h.GetNbinsX()+1):
		h.SetBinContent(i,j,i+10*j)
		print " i: ", i, " j: " ,j, " = ", i+10*j



		

h1 = unroll(h)

h1.SaveAs("h1.root")
h.SaveAs("h.root")

#	print "i: ", i
#	h.SetBinContent(i,-i)

'''
h.SetBinContent(0.5,0.5,-1)
h.SetBinContent(1.5,0.5,-2)
h.SetBinContent(1.5,0.f5,-3)
h.SetBinContent(1.5,0.5,-3)
h.SetBinContent(1.5,0.5,-3)
'''

'''


for y in range(0,h.GetNbinsY()+2):
	for x in range(0,h.GetNbinsX()+2):
		print " x : ", x , " y: ", y , " - ", h.GetBin(x,y), " - ", h.GetBinContent(h.GetBin(x,y))
'''

