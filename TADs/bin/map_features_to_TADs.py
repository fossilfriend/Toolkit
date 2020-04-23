#!/usr/bin/env python

"""
2019 Emily Greenfest-Allen

Description:
Generates index files mapping genomic features to TADs 

Features

1. Complete GENCODE gene reference
2. Variants (identified by chromosome position or chromosome start end)
3. Spans (e.g., custom gene list (with coordinates), LD windows)

----------
adapted from
2016 Gregory Way
scripts/generate_index_files.py
----------

Adaptations 
--> map through Ensembl IDS instead of gene symbols (names)
--> correct positional mismatch due to 1-based GTF files/variant coords and 0-based TAD bed files
--> improve readability w/meaningful variable names
--> improve efficiency of GTF file parsing
--> assign unmapped variants & spans to NaN instead of 0

"""

import argparse
import pandas as pd
import numpy as np
import gzip

from os import path
from Toolkit.Utils.utils import warning, qw, die

def parse_gzipped_gtf(gtf_file):
    '''
    Parameters:
    gtf_file: the full path to the gtf file (gzipped)
    Output:
    pandas dataframe w/gene info
    '''
    series = None
    lineCount = 0
    
    with gzip.open(gtf_file) as fh:
        for line in fh:
            if line.startswith('#'):
                continue
            lineCount = lineCount + 1
            line = line.rstrip()

            chromosome, db, featureType, start, stop, t1, strand, t2, info = line.split('\t')
            if featureType != 'gene':
                continue;
            gene_id, gene_name, gene_type = extract_relevant_gene_info(info)

            if series is None:
                series = [[chromosome.replace('chr', ''), int(start), int(stop), strand, gene_id, gene_type, gene_name]] # pd.Series([chromosome, start, stop, strand, gene_id, gene_type, gene_name])
            else:
                series.append([chromosome.replace('chr', ''), int(start), int(stop), strand, gene_id, gene_type, gene_name])
                
            if lineCount % 10000 == 0:
                warning("GENERATING GENE REF - Parsed", lineCount, "lines.")

        df = pd.DataFrame(series, columns = qw('chromosome start end strand gene_id gene_type gene_name', returnTuple=False))

        warning("GENERATING GENE REF - Done - Extracted", len(df.index), "genes.")
        
        return df
    

def extract_relevant_gene_info(infoStr):
    info = infoStr.split(';')
    build_info = {}
    for attribute in info:
        attribute = attribute.lstrip().replace('"', '')
        if attribute != '':
            attribute_name, new_attribute = attribute.split(' ')
            if attribute_name == 'gene_id':
                new_attribute = new_attribute.split('.')[0] # remove . in gene ids
            build_info[attribute_name] = new_attribute

    return [build_info['gene_id'], build_info['gene_name'], build_info['gene_type']]


def load_tads(TAD_file):
    '''
    Parameters:
    tad_file: full path to TAD file
    output: pandas dataframe containing the TADs
    '''
    warning("LOADING TADS")
    TAD_df = pd.read_csv(TAD_file, names = qw('chromosome start end'), sep='\t')
    TAD_df['chromosome'] = TAD_df['chromosome'].map(lambda x: x[3:])
    warning("LOADING TADS - Done")

    warning("TRANSFORMING TO 1-based to allow comparison to gene GTF and variant positional information")
    TAD_df['start'] = TAD_df['start'] + 1

    return TAD_df


def map_tad_elements(tads, features):
    '''
    Loop through TAD boundaries to assign genomic elements to TADs

    Arguments:
    :param tads: pandas dataframe, TAD boundary locations
    :param features: pandas dataframe, the genomic features to map to TADs

    Output:
    data frame with mapped tads to features

    '''

    result = pd.DataFrame()
    for tad in tads.itertuples():

        tad_id, chrom, start, end = list(tad)
        tad_id = 'tad_' + chrom + ':' + str(start) + '-' + str(end)

        # extract all overlapping spans
        overlappingFeatures = features[(features['chromosome'] == chrom) &
                                       (features['start'] <= end) &
                                       (features['end'] >= start)] # 1-based coordinate so should be inclusive
  
        if not overlappingFeatures.empty:
            overlappingFeatures['TAD'] = tad_id
            overlappingFeatures['TAD_start'] = start
            overlappingFeatures['TAD_end'] = end
     
            result = result.append(overlappingFeatures, ignore_index=True)

    return result



def flag_boundary_features(data):
    '''
    flag features that overlap TAD boundaries
    '''

    warning("FLAGGING TAD boundary crossing features")
    data['overlaps_TAD_boundary'] = (data['start'] < data['TAD_start']) | (data['end'] > data['TAD_end'])

    return data

                           
def map_genes(data, index_file):
    ''' 
    map genes to TAD regions
    '''
    
    warning("MAPPING Genes to TADs")
    result = map_tad_elements(tad_df, data)
    warning("MAPPING Genes to TADs - DONE")
       
    result = flag_boundary_features(result)
      
    result.to_csv(index_file, sep='\t', compression=args.compress, index=False)

   
def map_variants(inputFileName, index_file):
    '''
    map variants to TAD regions
    '''
    warning("MAPPING Variants to TADs")
    data = pd.read_csv(inputFileName, sep=args.delimiter)
    if 'position' not in data and 'start' not in data and 'end' not in data:
        die("Error parsing Variant data, expect either 'position', 'start',  or ['start', 'end'] columns")
    if 'position' in data:
        data['start'] = data['position']
        data['end'] = data['position']
    if 'start' in data and 'end' not in data:
        data['end'] = data['start']
    if 'chr' in data:
        data['chromosome'] = data['chr']
        
    data['chromosome'] = data['chromosome'].apply(str) # in case all were numeric (i.e., no X,Y)
    
    result = map_tad_elements(tad_df, data)
    warning("MAPPING Variants to TADs - DONE")
    
    unmapped = data.query("id not in @result.id")
    unmapped['TAD'] = np.nan
    unmapped['TAD_start'] = np.nan
    unmapped['TAD_start'] = np.nan

    result.append(unmapped, ignore_index=True, sort=False)
    
    result.to_csv(index_file, sep='\t', compression=args.compress, index=False)


def map_spans(inputFileName, index_file):
    ''' map spans to TAD regions '''
    warning("MAPPING TADs to Spans")
    data = pd.read_csv(inputFileName, sep=args.delimiter)

    if 'chr' in data:
        data['chromosome'] = data['chr']
    data['chromosome'] = data['chromosome'].apply(str) # in case all were numeric
    data['id'] = data['ld_block']
    
    result = map_tad_elements(tad_df, data)
    warning("MAPPING TADs to Spans - DONE")

    unmapped = data.query("id not in @result.id")
    unmapped['TAD'] = np.nan
    unmapped['TAD_start'] = np.nan
    unmapped['TAD_start'] = np.nan

    result.append(unmapped, ignore_index=True, sort=False)
    
    result.to_csv(index_file, sep='\t', compression=args.compress, index=False)


    
if __name__ == "__main__":
    '''generate index files'''

    pd.options.mode.chained_assignment = None

    parser = argparse.ArgumentParser()
    parser.add_argument('-c', '--cellType', help='boundary cell type', default='hESC')
    parser.add_argument('-t', '--tadFile', help='full path to 3-column-tab-separated TAD domain bed file', required=True)
    parser.add_argument('-i', '--inputFile', help="full path to input file (list of variants or spans); tab delim: must contain the following fields: 'id', 'chromosome', 'start', 'end'.  For genes, expects gzipped GENCODE GTF file.  For custom gene list, treat as spans", required=True)
    parser.add_argument('--featureType', help='gene, variant, or span', choices=['gene', 'variant', 'span'], default = 'variant', required=True)
    parser.add_argument('--prefix', help='prefix to add to output files')
    parser.add_argument('-o', '--outputFilePath', help='output file path', required=True)
    parser.add_argument('--genomeBuild', default='hg19')
    parser.add_argument('-d', '--delimiter', default='\t')
    parser.add_argument('--compress', help='compress output?', choices=['gzip', 'bz2'])
    args = parser.parse_args()

    # Read in TAD boundary file
    tad_df = load_tads(args.tadFile)

    extension = '' if not args.compress else \
        '.gz' if args.compress == 'gz' else '.bz2'
    
    if args.featureType == 'gene':
        # geneRef = generate_gene_ref(args.gencodeFile) if args.generateGeneIndex or args.featureType == 'ld' else None
        geneRef = parse_gzipped_gtf(args.inputFile) 
        indexFileName = args.prefix + '_gene_index_' + args.genomeBuild + '_' + args.cellType + '.tsv' + extension \
            if args.prefix else \
              'gene_index_' + args.genomeBuild + '_' + args.cellType + '.tsv' + extension 
        map_genes(geneRef, path.join(args.outputFilePath, indexFileName))
    
    if args.featureType == 'span':
        indexFileName = args.prefix + '_span_index_' + args.genomeBuild + '_' + args.cellType + '.tsv' + extension \
            if args.prefix else \
              'span_index_' + args.genomeBuild + '_' + args.cellType + '.tsv' + extension
     
        map_spans(args.inputFile, path.join(args.outputFilePath, indexFileName))
        
    if args.featureType == 'variant':
        indexFileName = args.prefix + '_variant_index_' + args.genomeBuild + '_' + args.cellType + '.tsv' + extension \
            if args.prefix else \
              'variant_index_' + args.genomeBuild + '_' + args.cellType + '.tsv' + extension
         
        map_variants(args.inputFile, path.join(args.outputFilePath, indexFileName))

