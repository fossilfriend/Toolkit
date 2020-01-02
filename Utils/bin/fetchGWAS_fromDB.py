#!/usr/bin/env python

'''
fetches GWAS data from GenomicsDB -- requires login credentials (gus.config) file

note -- only returns most significant p-value per variant

will update to query against webservice when permissions are monitored by saml

'''

from __future__ import print_function
from __future__ import with_statement

import argparse
import sys
from os import path
from Toolkit.Utils.postgres_dbi import Database
from Toolkit.Utils.utils import warning, xstr


SQL="""

WITH datasets AS (SELECT %s AS track)
SELECT DISTINCT
CASE WHEN split_part(metaseq_id, ':', 1) = 'X' THEN 23
WHEN split_part(metaseq_id, ':', 1) = 'Y' THEN 24
WHEN split_part(metaseq_id, ':', 1) = 'M' THEN 25
ELSE split_part(metaseq_id, ':', 1)::integer END AS chr,

split_part(metaseq_id, ':', 2)::integer AS pos,

allele AS testallele,

metaseq_id AS variant,
CASE WHEN r.source_id LIKE 'rs%%' THEN r.source_id  ELSE NULL END AS marker,

max(r.neg_log10_pvalue) OVER w AS neg_log10_pvalue,
first_value(r.pvalue_display::text) OVER w AS pvalue,
first_value(r.frequency) OVER w AS frequency,

first_value(r.restricted_stats->>'beta') OVER w AS beta,
first_value(power((r.restricted_stats->>'freq_se')::numeric, 2)) OVER w AS variance

FROM
Results.VariantGWAS r,
Study.ProtocolAppNode pan,
Datasets
WHERE pan.source_id = datasets.track
AND r.protocol_app_node_id = pan.protocol_app_node_id

WINDOW w as (PARTITION BY metaseq_id)


"""

def fetch_gwas_data():
    ''' fetch gwas data '''
    count = 0
    fileName = path.join(args.outputFilePath, args.track + ".txt")
    with open(fileName, 'w') as fh, database.cursor() as cursor:
        cursor.execute(SQL, (args.track,))
        for row in cursor:
            count = count + 1
            print('\t'.join((xstr(x) for x in row)), file=fh)
    warning("DONE - Wrote", count, "rows.")
                
if __name__ == "__main__":
    """fetch GWAS dataset"""
    parser = argparse.ArgumentParser("Query against NIAGADS GenomicsDB to extract GWAS result")
    parser.add_argument('-t', '--track', help="Unique track identifier", required=True)
    parser.add_argument('-o', '--outputFilePath', help="full output file path", required=True)
    parser.add_argument('--gusConfigFile',
                        help="GUS config file. If not provided, assumes default: $GUS_HOME/conf/gus.config")
    parser.add_argument('--verbose', action='store_true')

    args = parser.parse_args()

    
    database = Database(args.gusConfigFile)
    database.connect()

    fetch_gwas_data()
    
    database.close()
