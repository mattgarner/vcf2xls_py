"""Summary
"""
import os
import sys
import MySQLdb
import begin
import vcf
import xlwt
import urllib2

@begin.start
def decon_dict(decon_results_filepath):
    """Convert a DECoN results file into a python dict
    
    Args:
        decon_results_filepath (str): DECoN results fielpath
    
    Returns:
        dict: dict[cnv_id][field:value]
    """
    
    decon_data = {}
    cnv_id_header_text = "CNV.ID"

    with open(decon_results_filepath) as decon_fh:
        for line in decon_fh:
            line_components = line.strip().split("\t")
            
            if cnv_id_header_text in line_components:
                columns = line_components
                print columns
                continue

            cnv_id = int(line_components[columns.index(cnv_id_header_text)])
            
            # New entry for new CNVs
            decon_data.setdefault(cnv_id, {})

            kvps = {key:value for key, value in zip(columns, line_components)}

            for key, value in kvps.items():
                
                # Keep minimum start and max end values for each CNV
                if key == "Start.b":
                    value = int(value)
                    if value < decon_data[cnv_id].get(key, float("inf")):
                        decon_data[cnv_id][key] = value
                elif key == "End.b":
                    value = int(value)
                    if value > decon_data[cnv_id].get(key, -1):
                        decon_data[cnv_id][key] = value
                # Split multigenes and add any which are missing
                elif key == "Gene":
                    genes = value.split(";")
                    for gene in genes:
                        if gene not in decon_data[cnv_id].get("Gene", []):
                            decon_data[cnv_id].setdefault("Gene", []).append(gene)
                else:
                    # If no entry yet for this value in this cnv, add it
                    if decon_data[cnv_id].get(key, None) == None:
                        decon_data[cnv_id][key] = value
            

    for key in sorted(decon_data)[::-1]:
        print key
        for key, value in decon_data[key].items():
            print "%s: %s" % (key, value)
        print

    return decon_data