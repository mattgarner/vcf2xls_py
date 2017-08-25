import os
import sys
import MySQLdb
import begin
import vcf
import xlwt
import urllib2
import subprocess
import pprint
import re
import pysam
import pprint

# Data references
samplepanels_filepath = "/data/gemini/BioinformaticManifest.txt"
genepanels_filepath   = "/data/gemini/genepanels"
target_bed = "/data/projects/matt/vcf2xls_py/GENCODE_TSOE_v2_sorted.bed"

def sampleid2pids(sample_id, samplepanels_filepath):
    """Return a list of panel_ids and their genes for the specified g_number
    
    Args:
        g_number (str): Gemini sample identifier (e.g. G001234)
        samplepanels (TYPE): Filepath for tab delimited text file containing sample_id, Panel_Name, Panel_ID, Gene_Name
    
    Returns:
        dict: pid:panel_name for the specified sample
    """

    panel_ids = {}
    if sample_id.startswith("C"):
        pass

    elif any( [sample_id.startswith(prefix) for prefix in ["G","X","W"]]):
        with open(samplepanels_filepath) as samplepanels_fh:
            for line in samplepanels_fh:
                line_components = line.strip().split("\t")
                if len(line_components) != 4:
                    continue
                g_num, panel_name, pid, gene = line_components
                if g_num == sample_id:
                    # This is the data we need
                    # However, single gene panels have 'NA' as pid
                    # There should never be another version of these panels
                    # since they contain only one gene by definition,
                    # therefore for these we use panel_name as the identifier
                    if pid == "NA":
                        pid = panel_name

                    panel_ids[pid] = panel_name
        assert panel_ids, "No panels defined for %s in %s" % (sample_id, samplepanels_filepath)
        return panel_ids

def panelname2pid(panel_name, genepanels_filepath):
    """Return the panel_id for a single panel_name
    
    Args:
        panel_name (str): The name of a panel
        genepanels_filepath (str): Filepath to the panels data file
    
    Returns:
        str: A panel_id
    """
    if panel_name.startswith("_"):
        pid = panel_name
    else:
        with open(genepanels_filepath) as genepanels_fh:
            for line in genepanels_fh:
                line_components = line.strip().split("\t")
                
                # Skip incomplete records
                if len(line_components) != 3:
                    continue
                
                p_name, p_id, p_gene = line_components
                if panel_name == p_name:
                    pid = p_id
                    break
    return pid

def panelnames2pids(panel_names, genepanels_filepath):
    """Return a dict of panel_ids and panel_names
    
    Args:
        panel_names (list): A list of panel_names
        genepanels_filepath (str): Filepath to the panels data file
    
    Returns:
        dict: {pid:panel_name}
    """
    pid2panel = {}
    for panel_name in panel_names:
        pid = panelname2pid(panel_name, genepanels_filepath)
        pid2panel[pid] = panel_name
    return pid2panel

def get_sample_ids(sample_ids, vcf_in):
    vcf_sample_ids = list(vcf_in.header.samples)
    
    # If particular samples are specified we use those
    # otherwise use all samples in the vcf
    if sample_ids:
        sample_ids = sample_ids.replace(" ", "").strip().split(",")
        # All specified samples must be in the vcf
        for sample_id in sample_ids:
            assert sample_id in vcf_sample_ids, "Sample %s not found in %s" % (sample_id, vcf_filepath)    
    else:
        # All samples in vcf
        sample_ids = list(vcf_in.header.samples)

    return sample_ids

def sampleids2panels(sample_ids, panels):
    sample_to_panels = {}
    for sample_id in sample_ids:
        if panels:
            panel_names = panels.strip().split(",")
            panels = panelnames2pids(panel_names, genepanels_filepath)
            sample_to_panels[sample_id] = panels
        else:
            panels = sampleid2pids(sample_id, samplepanels_filepath)
            sample_to_panels[sample_id] = panels
    return sample_to_panels

def pid2genes(pid, genepanels_filepath):
    """Return a list genes for a specified panel_id
    
    Args:
        panel_ID (int): Unique identifier for panel(version specific)
        genepanels (str): Filepath for tab delimited text file containing Panel_Name Panel_ID Gene_Name
    
    Returns:
        list: Genes in the specified panel
    """

    panel_genes = []
    
    # Handle single gene panels with _GENE pids
    if pid.startswith("_"):
        gene = pid[1:]
        panel_genes.append(gene)

    # Handle full TSO
    elif pid == "8926":
        with open("/data/gemini/genes/gene_list_trusight_one.txt") as genelist_fh:
            for line in genelist_fh.readlines():
                gene = line.strip()
                panel_genes.append(gene)
    
    # Handle panels with numerical pids
    else:
        with open(genepanels_filepath) as genepanels_fh:
            for line in genepanels_fh:
                line_components = line.strip().split("\t")
                
                # Skip incomplete records
                if len(line_components) != 3:
                    continue
                
                panel_name, panel_ID, gene = line_components
                if pid == panel_ID:
                    panel_genes.append(gene)

    return list(set(panel_genes))

def pids2genes(pids, genepanels_filepath):
    """Return a dict of panel_ids and their genes
    
    Args:
        panels (list): A list of panel_IDs
        genepanels_filepath (str): Filepath for tab delimited text file containing Panel_Name Panel_ID Gene_Name
    
    Returns:
        dict: {pid:genes}
              genes contains an alphabetically sorted non-redundant list
              of names of all genes in all specified panels
    """
    
    # Some genes in the genepanels file have old HGNC identifiers. 
    # Consequently these gene names do not match those in the 
    # genes2transcripts file.
    # gene_alias_dict is used to map the old alias (the key) to the current
    # id in the genes2transcripts file(the value).

    panel_genes = {}

    #gene_alias_dict = get_gene_alias_dict()

    for pid in pids:
        genes = pid2genes(pid, genepanels_filepath)
        assert genes, "panel_id %s not found in genepanels_filepath" % pid
        # Replace old aliases
        #genes = [gene_alias_dict.get(gene, gene) for gene in genes]
        #missing_genes = get_missing_gene_list()
        #genes = [x for x in genes if x not in missing_genes]
        panel_genes[pid] = genes

    return panel_genes

def genes2regions(unique_genes, target_bed):
    # Gene -> Tx -> Exon -> Region
    with open(target_bed) as bed_fh:
        for line in bed_fh:
            chrom, start, end, region_id = line.strip().split("\t")
            gene_transcript_
    pass

@begin.start
def main(vcf_filepath, sample_ids=None, panels=None, maternal=None, paternal=None, index=None, relative=None):

    # Either single sample analyses, or a trio
    trio_analysis = all([maternal, paternal, index])
    individual_analysis = not any([maternal, paternal, index, relative])

    assert trio_analysis or individual_analysis, "\nFor trio analysis provide a vcf and maternal, paternal and index.\nFor individual analysis do not use maternal, paternal or index"

    # Check that all input files exist
    filepaths = [f for f in [vcf_filepath] if f]
    for f in filepaths:
        assert os.path.exists(f), "File not found at %s" % f

    vcf_in = pysam.VariantFile(vcf_filepath)
    
    sample_ids = get_sample_ids(sample_ids, vcf_in)
    
    # Only read the target samples - faster
    vcf_in.subset_samples(sample_ids)

    sample_panels = sampleids2panels(sample_ids, panels)

    # Sorry future me. This collapses the list of lists of pids. i.e. [[a,b,c],[a,d],[b,e]] -> [a,b,c,d,e]
    unique_pids = list(set([pid for pid_list in [x.keys() for x in sample_panels.values()] for pid in pid_list]))
    
    panel_genes = pids2genes(unique_pids, genepanels_filepath)

    unique_genes = list(set([gene for gene_list in panel_genes.values() for gene in gene_list]))

    gene_regions = genes2regions(unique_genes)

    print "Processing:"
    print "vcf:\t\t%s" % vcf_filepath
    pprint.pprint(sample_panels)
    pprint.pprint(panel_genes)

    #print dir(vcf_in)


'''
def main(
    vcf_filepath,
    decon_filepath="/data/projects/matt/cnv/DECoN_results/170510_170512_NS500192/170510_170512_NS500192_all.txt",
    samplepanels_filepath="/data/gemini/BioinformaticManifest.txt",
    genes2transcripts_filepath="/data/gemini/genes2transcripts",
    genepanels_filepath="/data/gemini/genepanels",
    panels=None,
    exon_flank=100,
    #text_output_dir="/data/projects/matt/vcf2xls_py/test_data/multialts/"  # DEBUGGING - REMOVE
    ):

    print "### Checking input ###"
    
    # Check that all input files exist
    input_filepaths = [ vcf_filepath,
                        samplepanels_filepath, 
                        genes2transcripts_filepath,
                        genepanels_filepath ]

    
    # Would be nice to auto detect this - requires a sensible storage system
    if decon_filepath:
        input_filepaths.append(decon_filepath)
    

    for filepath in input_filepaths:
        assert os.path.exists(filepath),\
            "File not found: %s" % filepath

    run_folder = os.path.basename(os.path.dirname(os.path.dirname(vcf_filepath)))
    sample_id = os.path.basename(vcf_filepath).split(".")[0]

    # Check that the vcf file contains data for the sample specified 
    # in the filename
    vcf_sample_ids = get_vcf_sample_ids(vcf_filepath)
    
    assert len(vcf_sample_ids) == 1,\
        "Multiple samples in vcf file"

    assert sample_id in vcf_sample_ids,\
        "SampleID in filename not found in file!"

    output_filepath = get_output_filepath(vcf_filepath)

    # Annotate any unannotated vcfs
    #vcf_filepath = annotate_vcf(vcf_filepath)

    #### Generate sample/panel specific data ####
    # NEED TO HANDLE CP numbers
    print "### Generating panel data ###"
    if panels:
        if panels == "full":
            panel_IDs = ["8926"]
        else:
            panel_IDs = panelnames2pids(panels, genepanels_filepath)
    else:
        panel_IDs = sampleid2pids(sample_id, samplepanels_filepath)

    gene_data_dict = get_gene_data_dict(sample_id, panel_IDs, genepanels_filepath, genes2transcripts_filepath, decon_filepath)

    panel_names = pids2panelnames(panel_IDs, genepanels_filepath)

    #### Generate vcf specific metadata ####
    csq_field_keys = get_csq_field_keys(vcf_filepath)

    # Set up the xls workbook
    with open(vcf_filepath) as vcf_fh:
        # strict_whitespace=True splits on tabs only, enabling spaces
        # to be present in fields which shouldn't really have them 
        # according to vcf spec.
        # Our vcfs do have such spaces so this is required to correctly
        # parse the vcfs
        vcf_reader = vcf.Reader(vcf_fh, strict_whitespace=True)

        styles = generate_styles()

        # Make a workbook
        wb = xlwt.Workbook()

        # Add the non-variant sheets to the workbook
        print "### Adding QC sheets ###"
        wb = add_summary_sheet(wb, sample_id, run_folder, styles, gene_data_dict, panel_names)
        wb = add_geneqc_sheet(wb, sample_id, styles, gene_data_dict)
        wb = add_exonqc_sheet(wb, sample_id, styles, gene_data_dict)

        print "### Setting up variant sheets ###"
        
        # Sheets and filters
        # For each variant sheet, provide a sheet name and a filtering
        # function which takes a row and returns true if the row should
        # be assigned to the sheet.
        # These functions will be executed in the order of the list, and
        # the row assigned to the sheet corresponding to the first 
        # function to return True. The final sheet function should
        # always return True to capture any variants not assigned to 
        # other sheets

        # TO DO - check each filter matches appropriate csq values in GO
        sheets = [  ("stop_gained", stop_gained_filter),
                    ("frameshift_variant", frameshift_filter),
                    ("splice", splice_filter),
                    ("missense", missense_filter),
                    ("synonymous", synonymous_filter),
                    ("other", other_filter),
                 ]

        sheet_names = [sheet[0] for sheet in sheets]

        # Specify column headers, functions to generate column values,
        # and col widths (in chars)
        # Cols will appear in the same order in the xls report
        # To add a new column add it to this list

        column_functions =     [("Gene", gene_col, 14),
                                ("Transcript", transcript_col, 16),
                                ("Chr", chromosome_col, 4),
                                ("Start", start_col, 11),
                                ("Ref", ref_col, 5),
                                ("Alt", alt_col, 5),
                                ("Nucleotide pos", c_pos_col, 12),
                                ("Nucleotide change", c_change_col, 12),
                                ("AA change", aa_change_col, 10),
                                ("Score", qual_score_col, 9),
                                ("Depth", depth_col, 7),
                                ("AAF", aaf_col, 6),
                                ("Genotype", genotype_col, 9),
                                ("dbsnp", dbsnp_col, 12),
                                ("Polyphen", polyphen_col, 14),
                                ("SIFT", sift_col, 14),
                                ("Gemini\nAF", af_gemini_col, 9),
                                ("Gemini\nHet/Hom", het_hom_gemini_col, 12),
                                ("ExAC\nAF", af_exac_col, 12),
                                ("ExAC\nHet/Hom", todo_col, 12),
                                ("Gnomad(E)\nAF", af_gnomad_exome_col, 12),
                                ("Gnomad(G)\nAF", af_gnomad_genome_col, 12),
                                ("Consequences", consequences_col, 30),
                                ("OMIM\nIM", OMIM_col, 10),
                                ]

        column_names = [x[0] for x in column_functions]
        column_widths = [x[2] for x in column_functions]
        #write_row_to_text_file(output_filepath, column_names)    # DEBUGGING - REMOVE

        # Add variant sheets to the workbook
        variant_sheets = {sheet[0]:wb.add_sheet(sheet[0]) for sheet in sheets}

        # Track the current row in each of the sheets to enable 
        # us to write to the correct line as we write records to 
        # different sheets
        current_row_index = {sheet[0]:0 for sheet in sheets}

        # Add variant sheet headers
        header_style = styles["header"]
        for sheet_name, ws in variant_sheets.items():

            # Each cell is written individually, so go through each col
            # in the row. For each cell we need a row index, a col index,
            # a value and a style
            for col_index, col in enumerate(column_functions):
                col_name, column_function, _ = col
                row_index = current_row_index[sheet_name]
                ws.write(row_index, col_index, col_name, style=header_style)

            freeze_header(ws)  # Keep header on screen when scrolling
            set_col_widths(ws, column_widths)
            header_row_height = 3
            set_row_height(ws, 0, header_row_height)
            current_row_index[sheet_name] += 1

        # Now populate the sheets with variants
        print "### Populating variant sheets ###"

        # Track the gene of the most recently added record in each sheet
        # Used to enable addition of empty lines between genes 
        current_gene_per_sheet = dict.fromkeys(sheet_names, None)

        for vcf_record in vcf_reader:

            # Get a list of alt alleles in the record
            alts = vcf_record.ALT

            # Get IDs of all transcripts of interest which the variant hits
            transcripts = record_in_regions(vcf_record, gene_data_dict, exon_flank)

            # Split csqs according to alt and transcript
            csq_dict = get_csq_dict(vcf_record, gene_data_dict, csq_field_keys)

            # Each transcript/alt combo generates a row
            for alt in alts:

                for transcript in transcripts:

                    row_data = []

                    # Get the relevant csq for the current 
                    # transcript/alt combo from the csq dict
                    csq = csq_dict[alt][transcript]

                    # Generate field values for the tx + alt combo
                    # by applying each col func to the data
                    for c in column_functions:
                        column_header = c[0]
                        column_function = c[1]
                        column_name, column_value = \
                            column_header,\
                            column_function(vcf_record=vcf_record,
                                            alt=alt, 
                                            transcript=transcript, 
                                            csq=csq,
                                            csq_field_keys=csq_field_keys)
                        row_data.append(column_value)

                    # Apply style rules
                    # Default:
                    row_style = styles["plain_text"]
                    # Exceptions
                    if float(row_data[column_names.index("Gemini\nAF")]) < 0.01:
                        row_style = styles["red_text"]
                    # Need purple rule

                    # Assign variant to a sheet using filter functions
                    for sheet in sheets:
                        sheet_name, sheet_filter = sheet
                        if sheet_filter(csq, csq_field_keys):
                            ws = variant_sheets[sheet_name]

                            # Add an empty row between genes
                            if (current_gene_per_sheet[ws.name] != None) \
                            and (row_data[column_names.index("Gene")] != current_gene_per_sheet[ws.name]):
                                current_row_index[ws.name] += 1 

                            write_row_to_worksheet(ws, current_row_index[ws.name], row_data, row_style)
                            #if len(alts) > 1:
                            #    print "\t".join(map(str,row_data))

                            current_row_index[ws.name] += 1
                            current_gene_per_sheet[ws.name] = row_data[column_names.index("Gene")]
                            break

    
    wb.save(output_filepath)                
    print "### Report complete! ###"
    print "Output: %s" % output_filepath
'''