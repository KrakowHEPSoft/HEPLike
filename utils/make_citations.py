#!/bin/python
import sys
import argparse
import os
import yaml



def main(argv):
    if(len(argv) != 2):
        print "Usage: python make_citations.py FILE.txt"
        print "where FILE.txt contains the bibcite names"

    f = open(argv[1], "r")
    lookup=[]
    for ff in f:
        lookup.append(ff.strip())
    fout = open("references.bib", "w")
    citations=[]
    #print lookup
    dir=""
    if(os.path.isdir("../data")):
        dir="../data"
    else:
        dir="data"
            
    for subdir, dirs, files in os.walk(dir):
        for file in files:
            if file.endswith(".yaml"):
                #yaml1 = yaml.YAML()
                #print 'Doing: ', subdir+"/"+file
                document = open(subdir+"/"+file)
                #print yaml.load(document)
                yamlf= yaml.load(document)
                #print yamlf['BibCite'] 
                if(yamlf['BibCite'] in lookup):
                    citations.append(yamlf['BibEntry'])
                    # we found already citation for this, we can remove it from list:
                    lookup.remove(yamlf['BibCite'])                    
                    # we found the citation, we can break
                    continue
                # now checking names:
                for l in lookup:
                    if( l in yamlf['Name']):
                        citations.append(yamlf['BibEntry'])
                        # we found the citation, we can break
                        continue
                    


    for c in citations:
        fout.write(c+'\n')
    fout.close()    

if __name__=="__main__":
    sys.exit(main(sys.argv))





