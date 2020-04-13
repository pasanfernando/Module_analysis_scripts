import networkx as nx
from collections import defaultdict
from operator import itemgetter
from collections import OrderedDict
import numpy as np

import curve_plotter_single_function as cp

#import matplotlib.pyplot as plt

""" This code is for module prediction phase
    predicts score for all the genes in the network
    Normalization is done during the code
    do not remove the particular gene in consideration from the gene list unlike the evaluation

    Code to gene function prediction based on the hishigaki methods
    This code is for evaluation using leave one out method
    The evaluation method for multiple classes/ functions is newly proposed
    
    Note: if you feel there is some thing wrong with the code, please uncomment the print statements throughout 
    the code. It will print results for each step. Do this with the prototype.
    """

############################################################################# ###########


#methods for calculating evaluation metrics

def evaluation_metrics(p,n,tp,tn,fp,fn):

    # calculating the accuracy
    try:
        accuracy = (tp+tn)/float(p+n)
    except ZeroDivisionError:
        accuracy = 0
    #print 'Accuracy is:',accuracy

    #error rate
    try:
        error_rate = (fp+fn)/float(p+n)
    except ZeroDivisionError:
        error_rate = 0
    #print 'Error rate is:', error_rate

    #sensitivity/ same as recall
    try:
        sensitivity = tp/float(p)
    except ZeroDivisionError:
        sensitivity = 0
    #print 'sensitivity:', sensitivity

    #specificity1
    try:
        specificity = tn/float(n)
    except ZeroDivisionError:
        specificity = 0
    #print 'specificity:', specificity

    #precision
    try:
        precision = tp/float(tp+fp)
    except ZeroDivisionError:
        precision=0
    #print 'precision:', precision

    #true positive rate
    try:
        tpr = tp/float(p)
    except ZeroDivisionError:
        tpr=0

    # false positive rate
    try:
        fpr = fp/float(n)
    except ZeroDivisionError:
        fpr=0

    # false positive rate

    return accuracy,error_rate,sensitivity,specificity,precision,tpr,fpr


########################################################################################

#input1 = raw_input('enter the interaction file:')
#input2 = raw_input('enter the genelist file:')

# prompting the user for gene list selection methody
input1 = raw_input('Extract a custom gene list from a file (y or n):')

if input1 == 'y':
    # a list to store the genes from the custom gene list
    clist =[]

    # if yes please check the name of the file
    cusgenelist = open('hindlimbgenelist_including_parts.txt', 'r')

    for line in cusgenelist:
        if line != '\n':
            line = line.strip()
            clist.append(line)

    cusgenelist.close()

# creating a nested dictionary to store chi square values for different functions for each gene
chigenes = defaultdict(dict)

#Reading the input gene function annotation file

genelist = open('newmouse_anatomy_profiles.txt', 'r')

# defining a multiple dictionary to store functions: genes
genefunc = defaultdict(list)
# defining a multiple dictionary to store gene: functions
genedic = defaultdict(list)

# Then convert the above input file into function: gene format

# this indicator is required for leave one out method

# reading the gene file
# for line in genelist:
#     if line != '\n':
#         line = line.strip()
#         a = line.split('\t') # splitting the line by tab to separate gene name vs functions
#         if a[0] != 'genes': # excluding the header
#             a[0] = a[0].lower() # convert to lower case to avoid mismatches
#             #print a[0]
#             # here, we are evaluating, so only the proteins with functions will be considered
#             if len(a)>1: # ths is required because there are proteins without one single function
#                 b = a[1].split(',') # seperating the functions by comma
#                 # store each function in the multi dic by function as the key and genes as values
#                 for i in b:
#                     genedic[a[0]].append(i) # this dictionary stores genes as keys
#                     genefunc[i].append(a[0])

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

# print genedic
# print len(genefunc)

### reading the network file########################################################
in1 = open('mslick0.3integrated_network.txt', 'r')
G = nx.Graph()
for line in in1:
    if line != '\n':
        #print line
        if 'combined_score' not in line:
            line = line.strip('\n')
            a = line.split()
            #print a[0],a[1],a[2]
            G.add_edge(a[0].lower(),a[1].lower(), weight = a[2])

# visualizing the network; please turn this off when loading large networks
# nx.draw_networkx(G, with_labels=True)
# plt.draw()
# plt.show()

######### calculate the function frequency or phi #########
# total proteins in the network
totps = nx.number_of_nodes(G)

# total proteins for the given function
# funps = len(genelist) # you are assuming all genes matches the ones in the network; need to check one by one

 # customized gene function list
customfunctions = ['uberon_0002103']
#customfunctions = ['uberon_0005724']

genes_in_network = nx.nodes(G)
print 'genes in the network:',genes_in_network

chvals =[]
# gene counter for runtime purposes
genecounter =0

genes_not_found = []
for gene in genes_in_network:
    print 'gene is:',gene
    genecounter +=1
    print genecounter

    #dictionaries to store total and function counts, chi square value
    total={}
    ccount={}
    chisq = {}

    #counting the genes with given function in immediate neighborhood of gene in consideration of leave one out method

    # going through the different functions in the function list
    for i in customfunctions:
    #for i in genefunc:
        #print 'function is:',i
        # reading the genes from a custom list if the option is given
        if input1 == 'y':
            genelist = clist
        else:
            genelist = genefunc[i]
        #print genelist

        # # making sure the gene in question is removed from the gene list because we are assuming it does not have any function
        # if gene in genelist:
        #     # removing the gene by list comprehension; otherwise it will mutate the original dictionary
        #     genelist = [x for x in genelist if x != gene]
        #print genelist
        #print genefunc[i]

        # making sure all the genes in the input data file are in the network

        notfound = set(genelist) - set(genes_in_network)  # genes that are not found in the network
        #print 'genes that are not found in the network:', notfound
        found = set(genelist) & set(genes_in_network)  # genes that are found in the network
        #print 'genes that are found in the network:', found


        funps = len(found)

        fracps = float(funps) / totps

        #print fracps

        # counting the genes with given function in immediate neighborhood of gene in consideration in leave one out method
        totalc =0 # initially, total count is zero
        ccountc=0 # initially, neighborhood count is zero
        if gene in G:
            nb = nx.all_neighbors(G, gene)
            #print i
            for j in nb:
                #print j
                totalc=totalc+1 # counting the total neighbors; could have done this using len(nb) as well
                if j in genelist: # if the neighbor is in the gene list that neghbor has the function
                    ccountc=ccountc+1

            # storing the above calculated values in designated dictionaries
            total[i]=totalc
            ccount[i]=ccountc

            #calculate expected frequency values
            ef = fracps*totalc
            # calculating the chi-square value
            try: # this error handling is required to avoid zero division when ef=0
                chi = ((ccountc-ef)**2)/ef
            except ZeroDivisionError:
                chi = 0
            chisq[i] = chi
            chvals.append(chi)
            chigenes[gene] = chisq
        else:
            genes_not_found.append(gene)

    #print chisq
    #print total
    #print ccount

print chigenes
# print chigenes['p1']
print chigenes

# getting the min max of chvals
minv = min(chvals)
maxv=max(chvals)

# ordering each dictionary based on the highest function

for i in chigenes:
    chigenes[i]= OrderedDict(sorted(chigenes[i].items(),key=lambda t:t[1], reverse=True))
    #print chigenes[i]

# # Normalizing the chi values
# for i in chigenes:
#
#     func_chi = chigenes[i]
#
#     #print lastelement[-1]
#     for j in func_chi:
#         #print j
#         #performing the normalization
#         func_chi[j]= (float(func_chi[j])-minv)/float(maxv-minv)
#
#     chigenes[i]=func_chi
#
# print chigenes

# a list to store chivalues; used to get the maximum chi value
chivalues = []
############################################################################################
# a dictionary to find the maximum chi value for each function

chimaxfunc = defaultdict(list)

# writing the predicted gene functions in an output file in sorted order
out = open('newfullhindlimbincludingparts_functions.txt', 'wb+')
out.write('genes\tpredicted_functions\n')

for i in chigenes:
    out.write('%s\t'%(i))
    func_chi = chigenes[i]
    lastelement = list(func_chi.keys())
    #print lastelement[-1]
    for j in func_chi:
        #print j
        # loading the chi values for each function to find the maximum chi value for each function
        chimaxfunc[j].append(func_chi[j])
        chivalues.append(float(func_chi[j]))# appending chi values
        if j == lastelement[-1]:# this avoids printing a comma after the last element
            out.write('%s:%s' % (j, func_chi[j]))
        else:
            out.write('%s:%s,'%(j,func_chi[j]))

    out.write('\n')
out.close()

# printing the chi max dictionary
print chimaxfunc
