from math import fabs, sqrt, log10, sin, cos, acos
from ROOT import TFile, TChain, TTree

def areEqual( x1, x2):
    if( x1 == 0. or x2 == 0.):
        if( fabs(x1-x2) > 1.e-5):
            return False
    else:
        if (x1 * x2 < 0.):
            return False
        else:
            if (fabs(log10(fabs(x1))-log10(fabs(x2))) > 1.e-5):
                return False
    return True

def dPhiFrom2P( Px1, Py1, Px2, Py2):
    prod = Px1*Px2 + Py1*Py2

    mod1 = sqrt(Px1*Px1+Py1*Py1)
    mod2 = sqrt(Px2*Px2+Py2*Py2)
    
    cosDPhi = prod/(mod1*mod2)
    
    return acos(cosDPhi)


def deltaR( Eta1, Phi1,	Eta2, Phi2):
    Px1 = cos(Phi1)
    Py1 = sin(Phi1)

    Px2 = cos(Phi2)
    Py2 = sin(Phi2)

    dPhi = dPhiFrom2P(Px1,Py1,Px2,Py2)
    dEta = Eta1 - Eta2

    return sqrt(dPhi*dPhi+dEta*dEta)

def load_entries(tree):
    l=[]
    for ievt in range(tree.GetEntries()):
        tree.GetEntry(ievt)
        l.append((ievt, int(tree.evt)))

    l.sort(key=lambda x: x[1])

    return l

class CompareSpring15:
    def __init__(self):
        self.ref=""
        self.test=""
        
        self.cref=None
        self.ctest=None

        self.lref=[]
        self.ltest=[]

        self.ref_missing=[]
        self.test_missing=[]

        self.bad_events=[]
        self.bad_leptons=[]
        self.bad_met=[]
        self.bad_mvamet=[]
        self.bad_jets=[]

        self.good_events=[]
    
    def load(self, ref, test):
        self.ref=ref
        self.test=test
                
        #self.cref=TChain(ref.split(":")[1])
        #self.ctest=TChain(test.split(":")[1])

        #self.cref.Add(ref.split(":")[0])
        #self.ctest.Add(test.split(":")[0])

	self.fref=TFile(ref.split(":")[0])
	self.ftest=TFile(test.split(":")[0])

	self.cref=TTree()
	self.fref.GetObject(ref.split(":")[1], self.cref)

	self.ctest=TTree()
        self.ftest.GetObject(test.split(":")[1], self.ctest)

        self.lref=load_entries( self.cref)
        self.ltest=load_entries( self.ctest)

        
    def Compare(self):
        print self.cref.GetEntries(), self.ctest.GetEntries()

        self.ref_min_pt_1 = 9999999999.;
        self.ref_min_pt_2 = 9999999999.;

        self.test_min_pt_1 = 9999999999.;
        self.test_min_pt_2 = 9999999999.;  
        
        iref=0
        itest=0

	self.cref.LoadBaskets(2000000000)
	self.ctest.LoadBaskets(2000000000)
        
        while(iref<len(self.lref)):
            ref_id=self.lref[iref][1]
            
            while(itest<len(self.ltest)):
                test_id=self.ltest[itest][1]
  
                if(test_id<ref_id):
                    self.ref_missing.append(self.ltest[itest])
                elif(test_id >ref_id):
                    self.test_missing.append(self.lref[iref])
                    break
                else:
                    break
                itest+=1
                
            if(test_id != ref_id):
                iref+=1
                continue

            if(self.cref.pt_1 < self.ref_min_pt_1):
                self.ref_min_pt_1 = self.cref.pt_1;
            if(self.cref.pt_2 < self.ref_min_pt_2):
                self.ref_min_pt_2 = self.cref.pt_2;
                
            if(self.ctest.pt_1 < self.test_min_pt_1):
                self.test_min_pt_1 = self.ctest.pt_1;
            if(self.ctest.pt_2 < self.test_min_pt_2):
                self.test_min_pt_2 = self.ctest.pt_2;

            isGood = True

            self.cref.GetEntry(self.lref[iref][0])
            self.ctest.GetEntry(self.ltest[itest][0])

            # compare leptons
            if( not areEqual(self.cref.pt_1, self.ctest.pt_1) or
                not areEqual(self.cref.iso_1, self.ctest.iso_1) or
                not areEqual(self.cref.pt_2, self.ctest.pt_2) or
                not areEqual(self.cref.iso_2, self.ctest.iso_2)
            ):

                isGood = False;

                print "Event", self.lref[iref][1]
                print "    lep1: pt=",self.cref.pt_1, self.ctest.pt_1,\
                        ", eta=",self.cref.eta_1,self.ctest.eta_1,\
                        ", phi=",self.cref.phi_1,self.ctest.phi_1,\
                        ", iso=",self.cref.iso_1,self.ctest.iso_1,\
                        ", q=",self.cref.q_1,self.ctest.q_1
                print "    lep2: pt=",self.cref.pt_2, self.ctest.pt_2,\
                        ", eta=",self.cref.eta_2,self.ctest.eta_2,\
                        ", phi=",self.cref.phi_2,self.ctest.phi_2,\
                        ", iso=",self.cref.iso_2,self.ctest.iso_2,\
                        ", q=",self.cref.q_2,self.ctest.q_2

                refDR=deltaR(self.cref.eta_1, self.cref.phi_1, self.cref.eta_2, self.cref.phi_2)
                testDR=deltaR(self.ctest.eta_1, self.ctest.phi_1, self.ctest.eta_2, self.ctest.phi_2)
                print "    deltaR(lep1, lep2)=",refDR,testDR

                self.bad_leptons.append((self.lref[iref],self.ltest[itest]))
                
            # compare met
            if( not areEqual(self.cref.met, self.ctest.met)):

                isGood = False;

                print "Event", self.lref[iref][1]
                print "    met: |met|=",self.cref.met, self.ctest.met

                self.bad_met.append((self.lref[iref],self.ltest[itest]))
                
            iref+=1
            itest+=1
            
    def Print(self):
        
        print len(self.ref_missing), "Events missing in ref tree"
        print len(self.test_missing), "Events missing in test tree"
        
        print "Events missing in test tree:"
        print [x[1] for x in self.test_missing]

        print
        print "Events missing in ref tree:"
        print [x[1] for x in self.ref_missing]

        print
        print len(self.bad_leptons), "Events with bad leptons"
        print [x[0][1] for x in self.bad_leptons]

        print
        print len(self.bad_met), "Events with bad pfmet"
        print [x[0][1] for x in self.bad_met]

        print self.ref_min_pt_1, self.test_min_pt_1
        print self.ref_min_pt_2, self.test_min_pt_2
