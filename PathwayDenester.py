import sys
import os
import io
import numpy as np
import pandas as pd
import math
import argparse
import warnings
from fractions import Fraction

parser = argparse.ArgumentParser(description="Look for pathways that are found only because of another more significant one:")

parser.add_argument('pathway_list', metavar='paths_address', type=str, nargs='?',
                    help='TSV result from pathways analysis like gprofiler, just columns \"term_id\" - Pathway ID as in the GMT file;  \"term_name\" - Pathway Name, aesthetic only.; \"intersection_size\"  (optional), \"term_size\"  (optional), and \"p_value\" or \"adjusted_p_value\" - From pathway enrichment, used for sorting pathways.; \"intersection\": List of genes that are differentially expressed in pathway, separated by \",\", must be quoted in case of CSV.') 

parser.add_argument('gmt_file', metavar='gmt_address', type=str, nargs='?',
                    help='GMT file with all pathways, will be used to find which genes are in each enriched pathway')

parser.add_argument('--selected_gene_list', metavar='selected_gene_list', type=str, nargs='?', default='',
                    help='Alternatively, if your pathway list doesn\'t show the genes that enrich each pathway, a separate list of genes can be inputed. Those will be comparaed to the gmt file to estimate which pathway they are associated with. This is LESS accurate as it might miss other parameters that might have been involved on the processed that created the pathway list. Only first column will be read, accepts the same nomenclaure used in the gmt file. Please remove genes that are not significant beforehand.') 

parser.add_argument('--output_address', metavar='output_address', type=str,  nargs='?', default='',
                    help='TSV file name where results will be saved. Default is pathway_list name  + \'_filtered.tsv\'')

parser.add_argument('--to_test_threshold', type=float, default='0',
                    help='Pathway \'B\' will be tested against all pathways \'A\' that are more significant AND that the ratio of B-DEGs also found in A to B-DEGs is at least [to_test_threshold].')

parser.add_argument('--pval_treshold', type=float, default='0.05',
                    help='P-value thresold to exclude a pathway. Since each pathway is treated independently, multiple testing corrections shouldn\'t be applied.')

parser.add_argument('--tranlator_gene_names', type=str, nargs='?', default='',
                    help='If you add a file that can translate the IDs used in gmt files to another name I can translate the genes list of each pathway. Argument are 3 comma separated strings; translator file address, colname of IDs, colname of translated names')

args = parser.parse_args()
print(args)
if (args.pathway_list == None and args.gmt_file == None):
    sys.exit('error; no input files, type "--help" to see the manual\n')
if (args.pathway_list == None):
    sys.exit('error; no input enrichment results, type "--help" to see the manual\n')
if (args.gmt_file == None):
    sys.exit('error; no pathway reference file, type "--help" to see the manual\n')

# Will transfer all arguments here for code clarity:
pathway_list = args.pathway_list
gmt_file =  args.gmt_file
output_address =  args.output_address
to_test_threshold =  float(args.to_test_threshold)   #default=0.00
pval_treshold = float(args.pval_treshold)  # to change result from "keep" to "exclude", default=0.05
tranlator_gene_names = args.tranlator_gene_names
selected_gene_list = args.selected_gene_list

version = '3.5'
#v2 I format the script neatly as a function
#v2.1 minor comments
#v2.4 adjusted identation, added a "filter" column to give more info on near misses
#v2.5 checks each line with the ones above intead of below.
#v2.6 sorts by density
#v2.7 sort by p-value again; tries to substiture density rule for reciprocity rule
#v2.8 accepts ties in the calculations, this means fully nested pathways are not excluded always anymore
#v2.8.1 calculates checks for ratio_of_unexpected_DEGs
#v2.9 deals with differences in p-value cutoffs, and splits the filter column into filter and reciprocal. Column order was rearranged
#v3 uses from fractions import Fraction
#v3.2 changes default to_test_threshold to 0. Fix strange case where there are p-value ties in the enriched input. Solve ties by density then number of degs
#v 3.3 accepts csv again
#v3.4 improves error messages and warnings
#3.5 added utf-8 capabilities in windows, also added optional --selected_gene_list argument

if output_address == '':
    output_address = pathway_list + '_filtered_' +version+ '.tsv'


def append_dict(dict, name, item):  # add to dictionary, if entry already exists, vector-append to it instead of replacing
	if name not in dict:
		dict[name] = []
	dict[name].append(item)

# This funtion is similar to panda's read_table(), but it returns a dictionary indexed by values from column dict_key. Synonyms are appended to a list.
# This function is larger as it performs some data cleaning. # I intend to reduce it in next versions.
def read_file(address, skiprows = 0, colorder = [], minlength = 2, dict_key = 0, coment_line_char = '#', sep = '\t'):
    dict_key_int = True
    if(str(type(dict_key)) == "<class 'int'>"):
        dict_key_int = True
    elif(str(type(dict_key)) == "<class 'list'>"):
        dict_key_int = False
    else:
        print('Error: type([dict_key]) must be \'list\' or \'int\', is: ' + str(type([dict_key])))
    colorder_int = True
    if(str(type(colorder)) == "<class 'int'>"):
        colorder_int = True
    elif(str(type(colorder)) == "<class 'list'>"):
        colorder_int = False
    else:
        print('Error: type([dict_key]) must be \'list\' or \'int\', is: ' + str(type([colorder])))

    end_of_line = 1  #will be used in an if several lines bellow.
    table = {}


    fil = io.open(address, 'r', encoding="utf-8")
    for line in fil:
        if(skiprows!= 0): #skip headers
            skiprows += -1
        else:
            if(end_of_line == 1):
                if(line[-2] == '\r'):
                    end_of_line = -2  # \r\n
                elif(line[-1] == '\n'):
                    end_of_line = -1  # \n
                else:
                    end_of_line = 0  #will happen if file is .gz, but the zcat might not be working
            if (minlength < len(line)):  # taking empty lines out
                if (line[0] != coment_line_char):  # taking commented lines out
                    splitted = line[:end_of_line].split(sep)  # taking out the '\n'
                    vector = []
                    if(colorder_int):
                        i = colorder
                        if(len(splitted) <= i):
                            print('Error on file ' + address + ' splitted len <= ' + str(i))
                            print(splitted)
                        if(dict_key_int):
                            append_dict(table, splitted[dict_key], splitted[i])
                        else:
                            append_dict(table, sep.join([splitted[j] for j in dict_key]), splitted[i])
                    else:
                        if(colorder != []):
                            for i in colorder:
                                if(len(splitted) <= i):
                                    print('Error on file ' + address + ' splitted len <= ' + str(i))
                                    print(splitted)
                                vector.append(splitted[i])
                            if(dict_key_int):
                                append_dict(table, splitted[dict_key], vector)
                            else:
                                append_dict(table, sep.join([splitted[j] for j in dict_key]), vector)
                        else:
                            if(dict_key_int):
                                for i in range(len(splitted)):
                                    #if(i != dict_key):  #comment this if you want the keys to be redundantly included in the vector
                                    vector.append(splitted[i])
                                append_dict(table, splitted[dict_key], vector)
                            else:
                                for i in range(len(splitted)):
                                    #if(i not in dict_key): #comment this if you want the keys to be redundantly included in the vector
                                    vector.append(splitted[i])
                                append_dict(table, sep.join([splitted[j] for j in dict_key]), vector)
    if(address[-3:] != '.gz'):
        fil.close()
    return (table)

def comb(n, k):
    if(k > n):
        return 0
    else:
        return(Fraction(math.factorial(n) , (math.factorial(k) * math.factorial((n - k)))))   #Fraction will remove the rounding errors that were building up very fast!

#mathematicaly derived function
#latex can be visualized in https://editor.codecogs.com/
#\sum_{k=r}^{min(n,m)}\frac{\binom{N-m}{n-k}*\binom{m}{k}}{\binom{N}{n}}
def comb_comb_comb(degs_in_test, degs_in_intersection, intersection_size, size_test):
    p_sum = 0
    denominator = comb(size_test,intersection_size)  #out of all possible combinations
    for desired_number_of_degs in range(degs_in_intersection, min(intersection_size, degs_in_test)+1):  #intersection_size is a sum of the probs of finding: all degs, all-1 deg and 1 nondeg, all-2 c degs and 2 nondegs, ...
        p_sum += Fraction(comb((size_test-degs_in_test),(intersection_size-desired_number_of_degs))*comb(degs_in_test,desired_number_of_degs), denominator) #get number of combinations of X degs and intersection-X nondegs
    return(float(p_sum))  #convert to float just at the very last moment

#lets load our gmt file
gmt_data = read_file(gmt_file, skiprows = 0, colorder = [], minlength = 2, dict_key = 0, coment_line_char = '#', sep = '\t')

#list all_known_genes:
all_known_genes = set()
for pathway in gmt_data:
    for gene in gmt_data[pathway][0][2:]:
        all_known_genes.add(gene)

if pathway_list.split('.')[-1].lower() == 'tsv':
    input_pathways = pd.read_csv(pathway_list, quotechar='\"', quoting = 0, sep = '\t')
elif pathway_list.split('.')[-1].lower() == 'csv':
    input_pathways = pd.read_csv(pathway_list, quotechar='\"', quoting = 0, sep = ',')
else:
    warnings.warn('PathwayDenester expects a tab separated file containing at leats the following columns:  \"term_id\", \"term_name\", and \"p_value\"')
    input_pathways = pd.read_csv(pathway_list, quotechar='\"', quoting = 0, sep = '\t')


###########
# adjust some nomenclature inconsistencies from different PEA sources
#make colun matching easier by converting them all to lowercase
input_pathways.columns = [colname.lower() for colname in input_pathways.columns]
#also replace - and space for _
input_pathways.columns = [colname.replace('-', '_') for colname in input_pathways.columns]
input_pathways.columns = [colname.replace(' ', '_') for colname in input_pathways.columns]

input_pathways = input_pathways.rename(columns={"intersections": "intersection"})  # won't be a problem if the name was already correct
input_pathways = input_pathways.rename(columns={"pathway_id": "term_id"})
input_pathways = input_pathways.rename(columns={"pathway_name": "term_name"})
input_pathways = input_pathways.rename(columns={"pvalue": "p_value"})
if 'p_value' not in input_pathways.columns:
    input_pathways = input_pathways.rename(columns={"adjusted_p_value": "p_value"})  #if there isn't a column called p_value, look for 'adjusted_p_value', they sort ~ the same way
    input_pathways = input_pathways.rename(columns={"adj_p_value": "p_value"})
    input_pathways = input_pathways.rename(columns={"q_value": "p_value"})


####### Sort pathways in input file; by p-value, by density in case of tie, by number of degs in case of tie
#kind='stable' preserves order in case of ties
#sort by number of degs
if 'intersection_size' in input_pathways.columns:
    input_pathways = input_pathways.sort_values(by='intersection_size', ascending=False, kind='stable')
    if 'term_size' in input_pathways.columns:
        input_pathways["ratio"] = (input_pathways["intersection_size"] / input_pathways["term_size"])
        input_pathways = input_pathways.sort_values(by='ratio', ascending=False, kind='stable')
input_pathways = input_pathways.sort_values(by='p_value', ascending=True, kind='stable')  #this one is required! Pathways will only be tested against those above them

        



rejected_pathways = input_pathways[~input_pathways.term_id.isin(gmt_data)] #because I need to know their composition to find the intersection sizes
if len(rejected_pathways) > 0:
    print('Warning, These pathways were not in the gmt file:\n' + '\n'.join([line for line in rejected_pathways.term_id]))
approved_pathways = input_pathways[input_pathways.term_id.isin(gmt_data)]
del input_pathways

spliting_string = ','  #hard coded  #how are the DEGs separated in referencre gmt

#####################
# Alternatively, if the pathway list doesn't show the genes that enrich each pathway, a separate list of genes can be inputed. Those will be comparaed to the gmt file to figure which pathway they are associated with. This is less accurate as it might miss other parameters that might have been involved on the processed that created the pathway list. Only fist column will be read.
deg_genes_set = set()
if ('intersection' not in approved_pathways.columns):
    if selected_gene_list != '':
        warnings.warn('No \'intersection\' column found pathway list, using '+selected_gene_list+' to estimate intersections...')
        gene_list_file = open(selected_gene_list, 'r')
        for line in gene_list_file:
            splitted = line[:-1].split('\t')
            splitted = splitted[0].split(',')  #split by either tab or comma, irrelevant as I will just read the 1st column
            deg_genes_set.add(splitted[0])
        gene_list_file.close()
        approved_pathways["intersection"] = '' #create new empty str() column
        #now use gmt to pupulate new column
        for pathway_id in approved_pathways.term_id:
            approved_pathways.loc[approved_pathways['term_id'] == pathway_id ,'intersection'] = spliting_string.join((deg_genes_set & set(gmt_data[pathway_id][0][2:])))
        #check for empty intersections:
        for i in approved_pathways.term_id:
            if (approved_pathways[approved_pathways['term_id'] == pathway_id]['intersection'] == '').item():  #.item() to extract bool from pandas object
                warnings.warn('No selected genes found in :'+pathway_id+' pathway will be excluded')
                approved_pathways = approved_pathways[approved_pathways['term_id'] != pathway_id]  #remove line
    else:
        sys.exit('error; \'intersection\' column not found on pathway list file, also no file was provided on --deg_genes_set.\nType "--help" to see the manual\n')


#@# In case we don't have a column listing the differentialy expressed genes (DEGs) in each pathway, but have a list of them separately
#@#pathways_have_independent_cutoffs = True


pathways_dictionaries = [{} for a in approved_pathways.iterrows()]
path_rank = 0
#load genes in pathway dictionary
for index, path in approved_pathways.iterrows():
    path_id = path.term_id
    pathaway_data = gmt_data[path_id][0]  #I need to know its composition to find the intersection size of the pathways
    path_slice = approved_pathways[approved_pathways.term_id == path_id].iloc[0]
    ratio_of_unexpected_DEGs = round(float(len(set(path_slice.intersection.split(spliting_string)) - set(pathaway_data[2:])))/len(path_slice.intersection.split(spliting_string)), 3)
    if (ratio_of_unexpected_DEGs > 0.1):   #hard coded; threshold to eliminate pathways where a high fraction of its DEGs not in the original gmt file
        print(path_id + ' has a high fraction of its DEGs not in the original gmt file: ' + str(ratio_of_unexpected_DEGs*100)+'%')
    else:
        approved_degs = [i for i in path_slice.intersection.split(spliting_string) if i in set(pathaway_data[2:])]
        pathways_dictionaries[path_rank] = {'id': path_id, 'name': path_slice.term_name, 'p-value': path_slice.p_value, 'all_genes': set(pathaway_data[2:]), 'deg_list': approved_degs, 'density' : float(len(approved_degs))/len(set(pathaway_data[2:]))}
        path_rank += 1

pathways_dictionaries = pathways_dictionaries[:path_rank]  # remove empty rows left by excluded pathways

#initialize new columns
for line in range(len(pathways_dictionaries)):
    pathways_dictionaries[line]['result'] = 1
    pathways_dictionaries[line]['filter'] = 'keep'
    pathways_dictionaries[line]['vs'] = 'itself'
    pathways_dictionaries[line]['vsName'] = ''
    pathways_dictionaries[line]['top10'] = pathways_dictionaries[line]['deg_list'][0:10]
    pathways_dictionaries[line]['degs'] = set(pathways_dictionaries[line]['deg_list'])
    pathways_dictionaries[line]['all_genes'] = set(pathways_dictionaries[line]['all_genes'])
    pathways_dictionaries[line]['reciprocal'] = 1


if tranlator_gene_names != '':
    print('translating gene IDs...')
    splitted = tranlator_gene_names.split(',')
    if len(splitted) != 3:
        sys.exit('error, tranlator_gene_names should be 3 commaseparated values, example:\n\tfile.tsv,gene_id,gene_name\nHow I reveived it :\n\t' + tranlator_gene_names + '\n')
    gene_translator_file = io.open(splitted[0], 'r', encoding="utf-8")
    header = gene_translator_file.readline()[:-1]
    header_split = header.split('\t')
    if (splitted[1] in header_split) and (splitted[2] in header_split):
        split_char = '\t'
    else:
        header_split = header.split(',')
        if (splitted[1] in header_split) and (splitted[2] in header_split):
            split_char = ','
        else:
            header_split = header.split(';')
            if (splitted[1] in header_split) and (splitted[2] in header_split):
                split_char = ';'
            else:
                 sys.exit('error, couldn\'t find columns \'' + splitted[1] + '\' and \'' + splitted[2]  + '\', I tried \'\\t\', \',\', and \';\' as field separators\n')
    col_from_position = header_split.index(splitted[1])
    col_to_position = header_split.index(splitted[2])
    gene_translator_file.close()
    gene_translator = read_file(splitted[0], skiprows = 1, colorder = col_to_position, minlength = 2, dict_key = col_from_position, coment_line_char = '', sep = split_char)
    for line in range(len(pathways_dictionaries)):
        new_top_10 = []
        for gene_id in pathways_dictionaries[line]['top10']:
            new_top_10.append(gene_translator[gene_id][0])
        pathways_dictionaries[line]['top10'] = new_top_10



#########
#if testing:
#start = time.time()
#fully_nested_pairs = [[], [], []]
#different_cutoffs = [[0], [0]]  # the zero starts here just for me to have a refernce of comparison in the resulting file.

#now lets iteratively do it for every pathway
for current_line in range(1,len(pathways_dictionaries)): #makes no sense to test 1st line
    degs_in_current = len(pathways_dictionaries[current_line]['degs']) 
    size_current = len(pathways_dictionaries[current_line]['all_genes'])
    for test_line in range(current_line):   #"test_line" is the "boss" line, of best p-value. "current_line" is the less significant pathway that may be booted. Index of test_line is always smaller than current_line.
        if pathways_dictionaries[test_line]['filter'] == 'keep':  #only test against those that weren't excluded yet
            degs_in_test = len(pathways_dictionaries[test_line]['degs'])
            #####################################
            # #Counting DEGs in the intersection between two pathways
            # ! When using ordered queries (sorted by p-value), two pathways can have different cutoffs. (added v2.9)
            # >> Regardless of the p-value cutoff used in the more significant pathway, we can only measure the independence of the blue pathway enrichment, if we use the same cutoff that was used in the enrichment we are trying to validate. 
            #degs_in_intersection_naive = len(pathways_dictionaries[test_line]['degs'] & pathways_dictionaries[current_line]['degs']) 
            #degs_in_intersection_test =  len(pathways_dictionaries[test_line]['degs'] & (pathways_dictionaries[test_line]['all_genes'] & pathways_dictionaries[current_line]['all_genes']))   #these two will differ in case the p-value threshold for the different pathwys is different ( ordred_querry = TURE )
            degs_in_intersection_current = len(pathways_dictionaries[current_line]['degs'] & (pathways_dictionaries[test_line]['all_genes'] & pathways_dictionaries[current_line]['all_genes']))  
            #####################################
            if (degs_in_intersection_current > to_test_threshold*min(degs_in_current, degs_in_test)):  #adjustable threshold #i will test if more than (ratio) their degs are from someone else
                intersection_size = len(pathways_dictionaries[current_line]['all_genes'] & pathways_dictionaries[test_line]['all_genes'])
                size_test = len(pathways_dictionaries[test_line]['all_genes'])
                current_result = comb_comb_comb(degs_in_current, degs_in_intersection_current, intersection_size, size_current)
                reverse_result = comb_comb_comb(degs_in_test, degs_in_intersection_current, intersection_size, size_test)
                if current_result < pval_treshold:
                    if reverse_result > pval_treshold:
                        ############################
                        # test the reciprocity rule:
                        #sometimes 2 pathways are mutually dependent, so the comb test would say that both can exclude the other. Let's only remove a pathway if the reverse wouldn't be true.
                        pathways_dictionaries[current_line]['result'] = current_result
                        pathways_dictionaries[current_line]['reciprocal'] = reverse_result
                        pathways_dictionaries[current_line]['filter'] = 'exclude'
                        pathways_dictionaries[current_line]['vs'] = pathways_dictionaries[test_line]['id']
                        pathways_dictionaries[current_line]['vsName'] = pathways_dictionaries[test_line]['name']
                        break
                if current_result < pathways_dictionaries[current_line]['result']:
                    ############################
                    #########even if it won't be filtered out, I will keep the best result I could find
                    if (pathways_dictionaries[current_line]['filter'] == 'keep'):
                        pathways_dictionaries[current_line]['vs'] = pathways_dictionaries[test_line]['id']
                        pathways_dictionaries[current_line]['vsName'] = pathways_dictionaries[test_line]['name']
                        pathways_dictionaries[current_line]['result'] = current_result
                        pathways_dictionaries[current_line]['reciprocal'] = reverse_result
                #########
                #if testing:
                #    print('line ' + str(current_line) + ': ' + str(time.time() - start))
                #    if degs_in_intersection_test != degs_in_intersection_current:
                #        different_cutoffs[0].append(current_line)
                #        different_cutoffs[1].append(test_line)
                #    if intersection_size ==  len(pathways_dictionaries[current_line]['all_genes']):  #if pathways are fully nested
                #        reverse_result = comb_comb_comb(degs_in_test, degs_in_intersection_current, intersection_size, size_test)
                #        fully_nested_pairs[0].append(pathways_dictionaries[test_line])
                #        fully_nested_pairs[1].append(pathways_dictionaries[current_line])
                #        fully_nested_pairs[2].append([test_line, size_test, degs_in_test, current_line, size_current, degs_in_current, intersection_size, degs_in_intersection_current, current_result, reverse_result, ((degs_in_current/size_current) > (degs_in_test/size_test))])



#############
#save results:
out_file = io.open(output_address, 'w', encoding="utf-8")
out_file.write('pathway id\tname\tpvalue\tDEG Density\tresult\treciprocal\tfiltered\tVersus\tVersusName\ttop10 genes\n')
#out_file.write('\t'.join([pathways_dictionaries[0]['id'], pathways_dictionaries[0]['name'], str(pathways_dictionaries[0]['p-value']), '2', 'best']) + '\n') #print line one
for line in range(0, len(pathways_dictionaries)):
    out_file.write('\t'.join([pathways_dictionaries[line]['id'], pathways_dictionaries[line]['name'], str(pathways_dictionaries[line]['p-value']),  str(round(pathways_dictionaries[line]['density'], 5)),  f"{pathways_dictionaries[line]['result']:.4g}", f"{pathways_dictionaries[line]['reciprocal']:.4g}" , str(pathways_dictionaries[line]['filter']), pathways_dictionaries[line]['vs'], pathways_dictionaries[line]['vsName'], ','.join(pathways_dictionaries[line]['top10'])]) + '\n')

out_file.close()

print('PathwayDenester done')


