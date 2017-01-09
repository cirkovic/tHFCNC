from ROOT import *
import sys

fs = []
#ts = []

ch = TChain('FlatTree/tree')

for a in sys.argv[2:]:
    fn = 'root://eoscms.cern.ch/'+a
    #fs.append(TFile('root://eoscms.cern.ch/'+a))
    #ts.append(fs[-1].Get('FlatTree/tree'))
    #ch.Add(fs[-1])
    ch.Add(fn)

print sys.argv[1], ch.GetEntries()

