#!/usr/bin/env python2

from contact_analyzer import PDBAnalyzer
from time import sleep



PDBAnalyzer=PDBAnalyzer()
#x=PDBAnalyzer.parse_structure("test.pdb")
#print x

#x,y,z=PDBAnalyzer._find_chains("test.pdb","A","C")
#print y

#dist=PDBAnalyzer._find_neighbors("test.pdb","A","C",10)
#print dist

#shortest=PDBAnalyzer.find_shortest_contacts("test.pdb","A","C",5)
#print shortest

#PDBAnalyzer.shortest_to_file("test.pdb","C","B",5,"shortest_distance.txt")

#PDBAnalyzer.contacts_to_file("test.pdb","C","B",5,"distances.txt")

#x,y,z=PDBAnalyzer._find_chains("./unittest/test_file.pdb","C","A")
#print x,y,z

#dist=PDBAnalyzer._find_neighbors("./unittest/test_file.pdb", "C", "A",5)
#print dist

#shortest=PDBAnalyzer.find_shortest_contacts("./unittest/test_file2.pdb", "C", "A",5)
#print shortest

#contacts=PDBAnalyzer.find_contacts("./unittest/test_file2.pdb", "C", "A",5)
#print contacts


PDBAnalyzer.shortest_to_file("test.pdb", "C", "A",5,"test.out")


#PDBAnalyzer.shortest_to_file("test.pdb","A","C",5,"output.txt")