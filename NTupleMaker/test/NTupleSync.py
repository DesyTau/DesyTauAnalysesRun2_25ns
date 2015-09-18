import sys
from DesyTauAnalyses.NTupleMaker.CompareSpring15 import areEqual, CompareSpring15
from DesyTauAnalyses.NTupleMaker.CompareSpring15 import CompareSpring15 as Compare

print areEqual(0.,0.0001)
print areEqual(3, 3)
print areEqual(3, 3.000001)
print areEqual(3, 3.00001)

print areEqual(0.00005, -0.0000000005)


print sys.argv[0], sys.argv[1], sys.argv[2]

c = Compare()
c.load(sys.argv[1], sys.argv[2])

c.Compare()

c.Print()
