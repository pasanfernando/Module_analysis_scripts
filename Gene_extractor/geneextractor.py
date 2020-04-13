import networkx as nx
from collections import defaultdict


'''
Extract set of genes for the given function from the network
Only the direct genes though
'''


function = 'uberon_0000151'


#Opening anatomical profile file
genelist = open('zebrafish_anatomy_profiles.txt', 'r')

# defining a multiple dictionary to store functions: genes
genefunc = defaultdict(list)
# defining a multiple dictionary to store gene: functions
genedic = defaultdict(list)

# reading the gene file in implementation
for line in genelist:
    if line != '\n':
        line = line.strip()
        a = line.split('\t') # splitting the line by tab to separate gene name vs functions
        if a[1] != 'gene_name': # excluding the header
            a[1] = a[1].lower() # convert to lower case to avoid mismatches
            #print a[0]
            # here, we are evaluating, so only the proteins with functions will be considered
            if len(a)>1: # ths is required because there are proteins without one single function
                b = a[3].split(',') # seperating the functions by comma
                # store each function in the multi dic by function as the key and genes as values
                for i in b:
                    genedic[a[1]].append(i) # this dictionary stores genes as keys
                    genefunc[i].append(a[1])

genelist = genefunc[function]

# a list to store selected genes
selectedgenes =[]

# Opening the original network

# storing genes in the network
netgenes =[]

# counter to count the interactions in the extracted network
intcounter=0

in1 = open('linstring0.33_network.txt', 'r')
# opening the preprocessed network output file
pre_net = open(function+'stringdirect_gene_extraction.txt', 'wb+')
# writing the header line
pre_net.write('protein1 protein2 combined_score\n')
for line in in1:
    if line != '\n':
        # print line
        if 'combined_score' not in line:
            templist = []
            line = line.strip('\n')
            a = line.split()
            # print a[0],a[1],a[2]
            gene1 = a[0].lower()
            # print gene1
            gene2 = a[1].lower()
            netgenes.append(gene1)
            netgenes.append(gene2)


            # checking wether both genes are in the newtork
            if gene1 in genelist and gene2 in genelist:
                pre_net.write('%s %s %s\n' % (gene1, gene2, a[2]))
                selectedgenes.append(gene1)
                selectedgenes.append(gene2)
                intcounter+=1


pre_net.close()
in1.close()

unselected = set(genelist)- set(selectedgenes)

#not in the network
notinnetwork = set(unselected)-set(netgenes)

alonegenes = set(unselected)&set(netgenes)

print 'total number of genes:',len(genelist)
print 'number of selected genes:',len(set(selectedgenes))
print 'number of unselected genes:',len(unselected)
print unselected
print 'number of genes not in the network',len(notinnetwork)
print notinnetwork
print 'number of genes in the network but alone:',len(alonegenes)
print alonegenes

#opening the statistics file
statfile = open('extractionstats.txt','wb+')
statfile.write('Total number of genes: %s\n'%(len(genelist)))
statfile.write('Number of selected genes: %s\n'%(len(set(selectedgenes))))
statfile.write('Number of interactions in the extracted network: %s\n'%(intcounter))
statfile.write('Total number of genes: %s\n'%(len(genelist)))
statfile.write('Number of unselected genes: %s\n'%(len(unselected)))
statfile.write('Number of genes not in the network: %s\n'%(len(notinnetwork)))
for i in unselected:
    statfile.write('%s\n'%(i))

statfile.write('\n')
statfile.write('Number of genes in the network but alone: %s\n'%(len(alonegenes)))
for i in alonegenes:
    statfile.write('%s\n'%(i))

statfile.write('\n')