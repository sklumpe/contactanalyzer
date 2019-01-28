#!/usr/bin/env python2

import Bio.PDB
import collections
import numpy as np
from interface import PDBAnalyzer as Base


class PDBAnalyzer(Base):

    def __init__(self):
        history=[]


    def _list_index_finder(self,lst,a):
        result=[]
        for i in range(0,len(lst)):
            if lst[i]==a:
                result.append(i)
        return(result)



    def parse_structure(self,pdb):
        parser=Bio.PDB.PDBParser()
        structures=parser.get_structure(str(pdb).replace(".pdb",""),pdb)
        structure=structures[0]
        return structure

    def _find_chains(self,pdb,chain1,chain2):
        structure=self.parse_structure(pdb)
        chain_1=[]
        chains_2={}
        chain_2=[]
        for chain in structure:
            if chain.get_full_id()[2]==chain1:
                chain_1=chain
            else:
                chains_2.update({chain.get_full_id()[2]:chain})
                chain_2.append(chain)
        return(chain_1,chains_2,chain_2)

    def _find_neighbors(self,pdb,chain1,chain2,cutoff=10):
        chain_1, chains_2, chain_2 = self._find_chains(pdb,chain1,chain2)
        structure = self.parse_structure(pdb)

        chain_2_atoms=Bio.PDB.Selection.unfold_entities(chain_2,"A")
        chain_1_residues=chain_1.get_residues()
        chain_1_atoms=Bio.PDB.Selection.unfold_entities(chain_1,"A")
        ns=Bio.PDB.NeighborSearch(chain_2_atoms)
        close_residue_dict={}
        ns_result=[]
        for i in chain_1:
            residue_atoms=i.get_list()
            close_residues=[]
            for j in residue_atoms:
                ns_result=ns.search(j.coord, cutoff, "R")
            if ns_result:
                close_residues.append(ns_result)
            flat_list=[item for sublist in close_residues for item in sublist]
            close_residue_dict.update({i.get_resname()+str(i.get_id()[1]):flat_list})
        return close_residue_dict, structure


    def find_contacts(self,pdb,chain1,chain2,cutoff=10):
        chain_1, chains_2, chain_2 = self._find_chains(pdb, chain1, chain2)
        close_residue_dict, structure = self._find_neighbors(pdb, chain1, chain2, cutoff)
        chain_1_residues=chain_1.get_residues()
        distance_list=[]
        for i in chain_1_residues:
            chain_1_atoms=i.get_list()
            resid=i.get_resname()+str(i.get_id()[1])
            for j in close_residue_dict[resid]:
                residue_atoms=j.get_list()
                for k in residue_atoms:
                    for h in chain_1_atoms:
                        distance_vector=k.coord - h.coord
                        distance=np.sqrt(distance_vector.dot(distance_vector))
                        if distance < cutoff:
                            distance_list.append([i,j,h,k,distance])
        return distance_list, close_residue_dict

    def contacts_to_file(self,pdb,chain1,chain2,cutoff=10,filename="contacts.out"):
        distance_list, close_residue_dict = self.find_contacts(pdb, chain1, chain2, cutoff)
        combined=[]
        entry_testlist = [(str(item[0].get_resname() + str(item[0].get_id()[1]))) for item in distance_list]
        i_testlist = [item[1] for item in distance_list]
        for n in range(0,len(i_testlist)):
            combined.append([entry_testlist[n],i_testlist[n]])
        with open(filename,"w") as output_file:
            for entry in close_residue_dict:
                for i in close_residue_dict[entry]:
                    index_list=self._list_index_finder(combined,[entry,i])
                    for j in index_list:
                        element_print=distance_list[j]
                        output_file.write("{0}	{1}	{2}	{3}	{4}	{5}	{6}\n".format(
                            entry,
                            chain1,
                            i.get_resname() + str(i.get_id()[1]),
                            chain2,
                            element_print[2].get_name(),
                            element_print[3].get_name(),
                            element_print[4]))




    def find_shortest_contacts(self, pdb, chain1,chain2,cutoff=10):
        chain_1, chains_2, chain_2 = self._find_chains(pdb, chain1, chain2)
        close_residue_dict, structure=self._find_neighbors(pdb,chain1,chain2,cutoff)
        chain_1_residues=chain_1.get_residues()
        distance_list=[]
        for i in chain_1_residues:
            chain_1_atoms=i.get_list()
            resid=i.get_resname()+str(i.get_id()[1])
            for j in close_residue_dict[resid]:
                residue_atoms=j.get_list()
                for k in residue_atoms:
                    for h in chain_1_atoms:
                        distance_vector=k.coord - h.coord
                        distance=np.sqrt(distance_vector.dot(distance_vector))
                        i_testlist = [item[0] for item in distance_list]
                        j_testlist = [item[1] for item in distance_list]
                        combined=[]
                        for n in range(0,len(i_testlist)):
                            combined.append([i_testlist[n],j_testlist[n]])
                        if [i,j] in combined:
                            position=combined.index([i,j])
                            if distance < distance_list[position][4]:
                                distance_list[position]=[i,j,h,k, distance]
                        else:
                            distance_list.append([i,j,h, k, distance])
        return distance_list, close_residue_dict

    def shortest_to_file(self,pdb,chain1,chain2,cutoff=10,filename="shortest_distance.out"):
        distance_list,close_residue_dict=self.find_shortest_contacts(pdb,chain1,chain2,cutoff)
        combined1=[]
        with open(filename,"w") as output_file:
            for entry in close_residue_dict:
                for i in close_residue_dict[entry]:
                    entry_testlist=[(str(item[0].get_resname()+str(item[0].get_id()[1]))) for item in distance_list]
                    i_testlist=[item[1] for item in distance_list]
                    for n in range(0,len(i_testlist)):
                        combined1.append([entry_testlist[n],i_testlist[n]])

                    if [entry,i] in combined1:
                        position=combined1.index([entry,i])
                        element_print=distance_list[position]
                        output_file.write("{0}	{1}	{2}	{3}	{4}	{5}	{6}\n".format(
												entry,
												chain1,
												i.get_resname()+str(i.get_id()[1]),
												chain2,
												element_print[2].get_name(),
												element_print[3].get_name(),
												element_print[4]))
                    else:
                        output_file.write("{0}	{1}	{2}	{3}	{4}	{5}	{6}\n".format(
                                                entry,
                                                chain1,
                                                i.get_resname()+str(i.get_id()[1]),
                                                chain2,
                                                "N/A",
                                                "N/A",
                                                "N/A"))



if __name__ == '__main__':
    import unittest
    import os
    class contact_analyzerTest(unittest.TestCase):
        def setUp(self):
            self.PDBAnalyzer = PDBAnalyzer()
        def test_shortest_to_file(self):
            self.PDBAnalyzer.shortest_to_file("./unittest/test_set_complex.pdb", "C", "A",10,filename="test_set.out")
            with open("./unittest/test_set_complex.out", "r") as ref:
                refLines = ref.readlines()
                refLines = [x[:-1] for x in refLines]
            with open("test_set.out", "r") as test:
                testLines = test.readlines()
                testLines = [x[:-10] for x in testLines]
            length=max((len(refLines),len(testLines)))
            self.assertEqual(set(refLines),set(testLines))


    unittest.main(verbosity=2)