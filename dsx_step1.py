#! /usr/bin/python

import dsx_config, re



genome = open("%s/%s" %(dsx_config.path_box,dsx_config.path_genome))
genomedict = dsx_config.open_genome(genome) #read genome file and save as dictionary

print genomedict["Scaffold1"]