import networkx as nx
from collections import defaultdict
import numpy as np

'''
Finding the best threshold for module detection
Using number of ortholog maps
this script is for zebrafish modules
'''

# the function that the prediction is done
function = 'uberon_0000151'

# the threshold value for the prediction
#threshold = 0

# defining dictionaries to store variables for threshold
# for orthologus genes
ogenes ={}
#for total number of genes in the module
togenemod ={}

# for the number of predicted genes
totpredicted ={}

# total number of unselected
totunselected ={}
# for number of genes not in the network
nonet ={}

# for isolated genes
isolated ={}



# loading the ortholog mapping file

# a multi dictionary to map zebrafish gene symbols: mouse symbols
zebradic = defaultdict(set)

# a multi dictionary to map mouse gene symbols: zebrafish symbols
mousedic = defaultdict(set)

# Opening the orthology mapping file downloaded from the zfin database

orthfile = open('mouse_orthos_2018.06.26.txt', 'r')

for line in orthfile:
    if line != '\n' and line.startswith('ZDB-GENE'):
        sep = line.strip().split('\t')

        # populating the zebrafish and mouse orthology dictionaries
        zebradic[sep[1].lower()].add(sep[3].lower())
        mousedic[sep[3].lower()].add(sep[1].lower())

# prompting the user for gene list selection method
input1 = raw_input('Extract a custom gene list from a file (y or n):')

if input1 == 'y':
    # a list to store the genes from the custom gene list
    clist = []

    # if yes please check the name of the file
    cusgenelist = open('genelist_including_parts.txt', 'r')

    for line in cusgenelist:
        if line != '\n':
            line = line.strip()
            clist.append(line)

    cusgenelist.close()
    # reading the genes from a custom list if the option is given
    orgenelist = clist
# Otherwise read from the original profile
else:
    orgenelist = genefunc[function]

# reading the mouse original gene list
# a list to store the genes from the custom gene list
otherlist = []

# if yes please check the name of the file
othergenelist = open('forelimbgenelist_including_parts.txt', 'r')

for line in othergenelist:
    if line != '\n':
        line = line.strip()
        otherlist.append(line)

#finding the threshold range from he precision file ####################################

# dictionary to store precision values: threshold
thresholddic ={}
prefile = open('uberon_0000151precision_hash.txt', 'r')
for line in prefile:
    if line != '\n':
        a=line.strip().split('\t')
        thresholddic[float(a[0])]=float(a[1])

itrange= np.linspace(0,30,100)
ftrange = np.linspace(30,max(thresholddic.keys()),100)
trange=  np.append(itrange,ftrange)

# opening the statistics table file
out= open('threshold_variation_stats.txt','wb+')
out.write('Threshold\tTotal number of modulegens\tNumber of predicted genes\tNumber of predicted genes confirmed as orthologs\tOrtholog percentage\tNumber of unselected genes\tNumber of genes not found in  the network\tNumber of isolated genes\n')

for threshold in trange:
    print 'threshold is:',threshold
    # creating a nested dictionary to store chi square values for different functions for each gene
    chigenes = defaultdict(dict)

    predgenes ={}
    # getting the predicted candidate genes
    # reading the chigenes dictionary from predicted functions file
    # a list to store chivalues; used to get the maximum chi value
    chivalues = []
    ############################################################################################
    # a dictionary to find the maximum chi value for each function

    # store all the genes and threshold score for the given function
    chimaxfunc = {}
    chidicfile = open('newfullpecincludingparts_functions.txt', 'r')
    for line in chidicfile:
        if line != '\n':
            line = line.strip()
            firstsplit = line.split('\t')
            if firstsplit[0] != 'genes':
                chigenes[firstsplit[0]]={}
                secondsplit= firstsplit[1].split(',')
                #print secondsplit
                for i in secondsplit:
                    #print i
                    thirdsplit = i.split(':')
                    chigenes[firstsplit[0]][thirdsplit[0]]=float(thirdsplit[1])
                    # if the predicted score is higher or equal to the threshold store that gene
                    if thirdsplit[0] == function:
                        if float(thirdsplit[1])>= threshold:
                            predgenes[(firstsplit[0])] = float(thirdsplit[1])
                        chimaxfunc[firstsplit[0]]=float(thirdsplit[1])
                    #chivalues.append(float(thirdsplit[1]))




    chidicfile.close()

    # geting the final gene list by adding original and final gene list
    genelist = set(orgenelist)|set(predgenes.keys())

    onlypredicted = set(predgenes.keys()) - set(orgenelist)
    # a list to store selected genes
    selectedgenes =[]
    print len(genelist)
    # Opening the original network

    # storing genes in the network
    netgenes =[]

    # counter to count the interactions in the extracted network
    intcounter=0

    in1 = open('newlinstring0.33preprocessed.txt', 'r')
    # opening the preprocessed network output file
    #pre_net = open(function+'_lin4gene_extraction.txt', 'wb+')
    # writing the header line
    #pre_net.write('protein1 protein2 combined_score\n')
    for line in in1:
        if line != '\n':
            # print line
            if 'combined_score' not in line:
                templist = []
                line = line.strip('\n')
                a = line.split()
                # print a[0],a[1],a[2]
                #print a
                gene1 = a[0].lower()
                # print gene1
                gene2 = a[1].lower()
                netgenes.append(gene1)
                netgenes.append(gene2)


                # checking wether both genes are in the newtork
                if gene1 in genelist and gene2 in genelist:
                    #pre_net.write('%s %s %s\n' % (gene1, gene2, a[2]))
                    selectedgenes.append(gene1)
                    selectedgenes.append(gene2)
                    intcounter+=1


    #pre_net.close()
    in1.close()

    # validating the orthologs
    print len(selectedgenes)


    commongenes =[]

    for gene in onlypredicted:
        if gene in otherlist:
            commongenes.append(gene)

        elif gene in zebradic:
            if zebradic[gene] in otherlist:
                commongenes.append(gene)



    print 'number of orthologous genes:',len(commongenes)

    ogenes[threshold]= len(commongenes)


    totpredicted[threshold]= len(onlypredicted)


    unselected = set(genelist)- set(selectedgenes)

    #not in the network
    notinnetwork = set(unselected)-set(netgenes)

    alonegenes = set(unselected)&set(netgenes)

    #print 'total number of original genes:',len(orgenelist)
    #print 'number of selected genes:',len(set(selectedgenes))
    togenemod[threshold]=len(set(selectedgenes))
    # print 'number of unselected genes:',len(unselected)
    # print unselected
    totunselected[threshold]=len(unselected)
    # print 'number of genes not in the network',len(notinnetwork)
    # print notinnetwork
    nonet[threshold]=len(notinnetwork)
    # print 'number of genes in the network but alone:',len(alonegenes)
    # print alonegenes
    isolated[threshold]=len(alonegenes)
    if len(onlypredicted)==0:
        orthlogratio=0.0
    else:
        orthlogratio= len(commongenes)/float(len(onlypredicted))

    out.write('%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n'%(threshold,len(set(selectedgenes)),len(onlypredicted),len(commongenes),orthlogratio,len(unselected),len(notinnetwork),len(alonegenes)))
