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
parser.add_argument("-o","--outputfile", help="Specify PDB output file", 
                required=False, type=str, default="extract.pdb")


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
obatoms = []
for atom in mol:
    coords.append(atom.coords)
    obatoms.append(atom.OBAtom)

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

    near_residues = {}
    for obatom in ob.OBResidueAtomIter(obres):
        atomcoord = numpy.array([obatom.x(), obatom.y(), obatom.z()])
        dist = numpy.sqrt(numpy.sum((npcoords - atomcoord)**2, axis=1))
        print "    %10.3f %10.3f %10.3f %s"%(obatom.x(), \
                obatom.y(), obatom.z(),  obres.GetAtomID(obatom))

        for index in numpy.flatnonzero(dist < args.mindist):
            res = obatoms[index].GetResidue()
            uniqname = res.GetName()+ "_" + str(res.GetIdx()) + "_" + \
                    res.GetChain()
            near_residues[uniqname] = res

    print "Now copying original molecule and removing atoms"
    extract_pdb = ob.OBMol(obmol)
    atomtoremove = []
    for res in ob.OBResidueIter(extract_pdb):
        uniqname = res.GetName()+ "_" + str(res.GetIdx()) + "_" + \
                    res.GetChain()
        if not ( uniqname in near_residues):
            for obatom in ob.OBResidueAtomIter(res):
                atomtoremove.append(obatom)

    for atomtorm in atomtoremove:
        extract_pdb.DeleteAtom(atomtorm)

    extract_pdb.AddHydrogens()
        
    output = pybel.Outputfile("pdb", args.outputfile, overwrite=True)
    output.write(pybel.Molecule(extract_pdb))
    output.close()

else:
    print "There are ", numofresidue, " with the same name "
