#!/bin/python
import sys
import argparse
import os
import yaml
import argparse


def main(argv):
    parser = argparse.ArgumentParser()
    parser.add_argument("--HLAuthor", help="search for Author of HepLike measurement file")
    parser.add_argument("--Year", help="search for measurments in the following year")
    #parser.add_argument("--Particles", help="search for particles in decay descriptor")  
    parser.add_argument("--Name", help="search for name of measurements") 
    parser.add_argument("--Arxiv", help="search for arxiv number")     

    args = parser.parse_args()
    
    
    dir=""
    
    if(os.path.isdir("../data")):
        dir="../data"
    else:
        dir="data"

    files_found=[]
        
    for subdir, dirs, files in os.walk(dir):
        for file in files:
            if file.endswith(".yaml"):
                document = open(subdir+"/"+file)
                yamlf= yaml.load(document)
                if('Arxiv' in yamlf):
                    if(yamlf['Arxiv'] == args.Arxiv):
                        files_found.append(subdir+"/"+file)
                        continue
                if('HLAuthor' in yamlf):
                    if(yamlf['HLAuthor'] == args.HLAuthor):
                        files_found.append(subdir+"/"+file)
                        continue
                    elif(str(args.HLAuthor) in yamlf['HLAuthor']):
                        files_found.append(subdir+"/"+file)
                        continue
                    
                        
                if( ('SubmissionYear' in yamlf)  and ('PublicationYear' in yamlf)):
                    if(  (yamlf['SubmissionYear'] == args.Year)   or (yamlf['PublicationYear']  == args.Year)    ):
                        files_found.append(subdir+"/"+file)
                        continue
    print 'Found files:'
    for f in files_found:
        print f


if __name__=="__main__":
        sys.exit(main(sys.argv))    
    
