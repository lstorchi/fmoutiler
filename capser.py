import argparse
import pybel
import numpy
import sys

import os.path

import openbabel as ob
from scipy.spatial import distance

parser = argparse.ArgumentParser()
parser.add_argument("-i","--inputfile", help="Specify PDB input file", 
                required=True, type=str, default="test.pdb")
parser.add_argument("-l","--referenceligand", help="Specify the ligand within the protein", 
                required=True, type=str, default="NON")
parser.add_argument("-m","--mindist", help="Specify minimum distance from ligand", 
                required=False, type=numpy.float64, default=4.0)


args = parser.parse_args()

print "Options: "
print args 
print ""

if not os.path.isfile(args.inputfile):
    print "PDB File ", args.inputfile, " does not exist"
    exit(1)

ligandname = args.referenceligand

mol = pybel.readfile("pdb", args.inputfile).next()
coords = []
for atom in mol:
    coords.append(atom.coords)


obmol = mol.OBMol 

obres = ""
numofresidue = 0
for res in ob.OBResidueIter(obmol): 
    if res.GetName() == ligandname:
        numofresidue += 1
        obres = res

if numofresidue == 1:
    npcoords = numpy.array(coords)
    print "Found ", obres.GetName(), " ", obres.GetIdx(), \
            " chain: ", obres.GetChain()

    for obatom in ob.OBResidueAtomIter(obres):
        atomcoord = numpy.array([obatom.x(), obatom.y(), obatom.z()])
        dist = numpy.sqrt(numpy.sum((npcoords - atomcoord)**2, axis=1))
        print "    %10.3f %10.3f %10.3f %s"%(obatom.x(), \
                obatom.y(), obatom.z(),  obres.GetAtomID(obatom))

        print "    near atom ",numpy.flatnonzero(dist < args.mindist)

else:
    print "There are ", numofresidue, " with the same name "
