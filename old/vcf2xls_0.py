"""Summary
"""
import os
import sys
import MySQLdb
import begin
import vcf
import xlwt


def geminiDB_query(query):
    db = MySQLdb.connect(host="mgsrv01",
                         user="easih_ro",
                         db="GeminiDB")

    #  you must create a Cursor object. It will let
    #  you execute all the queries you need
    cur = db.cursor()

    # SQL
    cur.execute(query)

    result = []
    for row in cur.fetchall():
        result.append(row)

    db.close()
    return result


def gnum2pids(g_number, samplepanels_filepath):
    """Return a list of panel_ids and their genes for the specified g_number
    
    Args:
        g_number (str): Gemini sample identifier (e.g. G001234)
        samplepanels (TYPE): Filepath for tab delimited text file containing G_Number, Panel_Name, Panel_ID, Gene_Name
    
    Returns:
        dict: gnum2panels_dict[panel_id]{"name":panel_name, genes:[gene1, gene2...]}
    """
    
    panel_ids = []
    with open(samplepanels_filepath) as samplepanels_fh:
        for line in samplepanels_fh:
            line_components = line.strip().split("\t")
            if len(line_components) != 4:
                continue
            g_num, panel_name, pid, gene = line_components
            if g_num == g_number:
                # This is the data we need - populate the output dict
                
                # Some single gene panels have 'NA' as pid
                # There should never be another version of these panels
                # since they contain only one gene by definition. 
                # therefore for these so use panel_name as key instead
                if pid == "NA":
                    pid = panel_name

                panel_ids.append(pid)

    return list(set(panel_ids))


def pids2panelnames(pids, genepanels_filepath):
    panel_names = {}
    
    for pid in pids:
        panel_name = pid2panelname(pid, genepanels_filepath)
        panel_names[pid] = panel_name

    return panel_names


def pid2panelname(pid, genepanels_filepath):
    if pid.startswith("_"):
            panel_name = pid
    
    else:
        with open(genepanels_filepath) as genepanels_fh:
            for line in genepanels_fh:
                line_components = line.strip().split("\t")
                
                # Skip incomplete records
                if len(line_components) != 3:
                    continue
                
                panel_name, panel_id, gene = line_components
                if pid == panel_id:
                    panel_name = panel_name
                    break

    return panel_name


def pids2genes(panel_IDs, genepanels_filepath):
    """Return a dict of panel_ids and their genes
    
    Args:
        panels (list): A list of panel_IDs
        genepanels_filepath (str): Filepath for tab delimited text file containing Panel_Name Panel_ID Gene_Name
    
    Returns:
        TYPE: A alphabetically sorted non-redundant list of names of all genes in all specified panels
    """
    
    panel_genes = {}

    for pid in panel_IDs:
        genes = pid2genes(pid, genepanels_filepath)
        panel_genes[pid] = genes
    
    return panel_genes


def pid2genes(panel_ID, genepanels_filepath):
    """Return a list genes for a specified panel_id
    
    Args:
        panel_ID (int): Unique identifier for panel(version specific)
        genepanels (str): Filepath for tab delimited text file containing Panel_Name Panel_ID Gene_Name
    
    Returns:
        list: Genes in the specified panel
    """

    panel_genes = []
    
    # Handle single gene panels with _GENE pids
    if panel_ID.startswith("_"):
        gene = panel_ID[1:]
        panel_genes.append(gene)
    
    # Handle panels with numerical pids
    else:
        with open(genepanels_filepath) as genepanels_fh:
            for line in genepanels_fh:
                line_components = line.strip().split("\t")
                
                # Skip incomplete records
                if len(line_components) != 3:
                    continue
                
                panel_name, pid, gene = line_components
                if pid == panel_ID:
                    panel_genes.append(gene)
    
    return panel_genes


def genes2transcripts(genes, genes2transcripts_filepath):
    """Return a dict of gene_name:[transcripts] for all genes
    
    Args:
        genes (list): A list of HGNC gene names
        genes2transcripts_filepath (str): Filepath for tab delimited text file containing Gene_Name Transcript
    
    Returns:
        dict: gene_name:[transcripts]
    """

    gene_transcripts = {}
    
    for gene_list in genes.values():
        for gene in gene_list:
            transcripts = gene2transcripts(gene, genes2transcripts_filepath)
            gene_transcripts[gene] = transcripts
        
    return gene_transcripts


def gene2transcripts(gene, genes2transcripts_filepath):
    """Return a list of transcripts for a given gene
    
    Args:
        gene (str): HGNC gene name
        genes2transcripts (str): Filepath for tab delimited text file containing Gene_Name Transcript_ID
    
    Returns:
        list: Transcripts used of the gene
    """
    transcripts_for_gene = []
    
    with open(genes2transcripts_filepath) as g2t_fh:
        for line in g2t_fh:
            line_components = line.strip().split("\t")
            if len(line_components) != 2:
                continue
            gene_name, transcripts = line_components
            if gene == gene_name:
                # Some genes have two NM_s in one record in the database
                # Third normal form grumble grumble...
                transcripts = transcripts.split(",")
                for transcript in transcripts:
                    transcripts_for_gene.append(transcript)
    
    return transcripts_for_gene


def transcripts2exons(transcripts):
    transcript_exons = {}
    
    for gene, transcripts in transcripts.items():
        for transcript in transcripts:
            exons = transcript2exons(transcript)
            transcript_exons[transcript] = exons
    return transcript_exons


def transcript2exons(transcript):
    # TO DO +/- range on exons
    """Summary
    
    Args:
        transcript (TYPE): Description
    
    Returns:
        TYPE: Description
    """
    query = """SELECT r2g.exon_nr, reg.chr, reg.start, reg.end
                 FROM gene g, region2gene r2g, region reg 
                WHERE g.gid = r2g.gid 
                  AND r2g.rid = reg.rid 
                  AND refseq like '%%%s%%'""" % transcript
    
    data = geminiDB_query(query)
    exons = {}

    for row in data:
        exon_nr, chrom, start, end = row
        exons[str(exon_nr)] = {"chr":str(chrom),
                               "start":int(start),
                               "end":int(end)}
    return exons


def record_in_regions(vcf_record, regions):
    vcf_record_range = record_range(vcf_record)
    
    
    for transcript, exons in regions.items():
        for exon, coords in sorted(exons.items()):

            # Wrong chrom
            if vcf_record.CHROM != coords["chr"]:
                break
            
            # Right chrom - do coords overlap?
            else:
                region_range = [coords["start"], coords["end"]]
                
                overlap_length = overlap(region_range, vcf_record_range)
                if overlap_length:
                    print "v", vcf_record_range
                    print "p", region_range
                    print overlap_length
                    return True
    return False
                   

def record_range(vcf_record):
    alleles = vcf_record.ALT + [vcf_record.REF]
    max_allele_length = max([len(allele) for allele in alleles])
    
    first = vcf_record.POS
    last = vcf_record.POS + max_allele_length

    return [first, last]


def overlap(a, b):
    return max(0, min(a[1], b[1]) - max(a[0], b[0]))


def xls_generator():
    pass


def gene_QC_sheet():
    pass


def exon_QC_sheet():
    pass


def prepare_variant_sheet():
    pass


def populate_variant_sheet():
    pass

'''
OBSOLETE - SEE get_csq_field_keys()
def get_csq_index(vcf_filepath):
    with open(vcf_filepath) as vcf_fh:
        vcf_reader = vcf.Reader(vcf_fh)
    
        format_tag = "Format: "
        csq_desc = vcf_reader.infos.get("CSQ").desc
        csq_format_pos = csq_desc.find(format_tag)
        csq_fields = csq_desc[csq_format_pos+len(format_tag):].split("|")
        csq_fields_index = {}
        for index, field in enumerate(csq_fields):
            csq_fields_index[field] = index
        return csq_fields_index
'''

def get_csq_field_keys(vcf_filepath):
    """Return an ordered list of csq field keys
    
    Args:
        vcf_filepath (str): Filepath for vcf file to be analysed
    
    Returns:
        List: Consequence field names in the order in which they appear in the vcf
    """

    with open(vcf_filepath) as vcf_fh:
        vcf_reader = vcf.Reader(vcf_fh)
        csq_desc = vcf_reader.infos.get("CSQ").desc
        csq_field_keys = csq_desc.split("Format: ")[-1].split("|")
        return csq_field_keys


#### Report columns ####
# Any higher order functions need to be called later to generate their
# function after the necessary sample specific data has been generated

def generate_gene_col(gene_transcripts):
    def gene_col(vcf_record):
        csqs = vcf_record.INFO["CSQ"]
        for csq in csqs:
            gene=csq.split("|")[2]
            if gene in gene_transcripts.keys():
                return gene
    return gene_col


def transcript_col(vcf_record):
    pass

def chromosome_col(vcf_record):
    return vcf_record.CHROM


def start_col(vcf_record):
    # Does this need to be POS+1 for non-SNPs?
    return vcf_record.POS

def ref_col(vcf_record):
    return vcf_record.REF

def alt_col(vcf_record):
    # Need to split multialts?
    return vcf_record.ALT

def pos_col(vcf_record):
    # Does this need to be POS+1 for non-SNPs?
    return vcf_record.POS

def generate_nuc_change_col(csq_fields_index, gene_transcripts):
    def nuc_change_col(vcf_record):

        # Where in the csq string is the nuc_change field?
        csq_field = get_csq_field("HGVSc", vcf_record, gene_transcripts, csq_fields_index)
        nuc_change = csq_field.split(":")[1]
        return nuc_change
            
    return nuc_change_col

def generate_aa_change_col(csq_fields_index, gene_transcripts):
    def aa_change_col(vcf_record):
        # NEEDS AN AA POS TOO

        csq_field = get_csq_field("Amino_acids", vcf_record, gene_transcripts, csq_fields_index)
        aa_change = csq_field.split("/")
        return aa_change
            
    return aa_change_col

def qual_score_col(vcf_record):
    return vcf_record.QUAL

def depth_col(vcf_record):
    pass

def aaf_col(vcf_record):
    for tag in ["ABHet", "ABHom"]:
        aaf = vcf_record.INFO.get(tag)
        if aaf:
            return aaf

def genotype_col(vcf_record):
    return [call.data.GT for call in vcf_record.samples]


def dbsnp_col(vcf_record):
    # Need to check which dbSNP is used
    return vcf_record.ID

def polyphen_col(vcf_record):
    pass

def sift_col(vcf_record):
    pass

def af_gemini_col(vcf_record):
    # NEEDS CHANGING TO GENERATE DYNAMICALLY - THIS TABLE IS NOT UP TO DATE
    af_gemini = []
    for alt in vcf_record.ALT:
        query = """SELECT freq, homs, hets from variant_freqs
                   WHERE vid IN 
                   (SELECT vid 
                    FROM variant 
                    WHERE chr = '%s' 
                    AND pos = '%s' 
                    AND ref = '%s' 
                    AND alt = '%s' 
                    AND prid = 1)""" % (vcf_record.CHROM, vcf_record.POS, vcf_record.REF, alt)
        freq = geminiDB_query(query)
        af_gemini.append(freq)
    return af_gemini
    # Need to know:
    #   how many samples have coverage x
    #   how many samples have the variant as het y
    #   how many samples have the variant as het z
    #   AF = (y+2z)/2x
    #   hets - count or %, if % as % of total or % of variant carriers
    #   homs - as hets
    pass

def af_1kg_max_col(vcf_record):
    pass

def af_esp_max_col(vcf_record):
    pass

def af_exac_col(vcf_record):
    pass

def af_gnomad_exome_col(vcf_record):
    pass

def af_gnomad_genome_col(vcf_record):
    pass

def comment_col(vcf_record):
    pass

def OMIM_IM(vcf_record):
    # Use OMIM API to pull inheritance mode info for gene
    pass


def get_csq_field(csq_field, vcf_record, gene_transcripts, csq_field_keys):
    # Needs to handle multi alts
    target_field_index = csq_field_keys.index(csq_field)
    csqs = vcf_record.INFO.get("CSQ")

    target_transcript_csqs = []
    for csq in csqs:
        transcript_field_name = "RefSeq"
        transcript_field_index = csq_field_keys.index(transcript_field_name)
        
        csq_fields = csq.split("|")
        csq_transcripts = csq_fields[transcript_field_index]
        
        for gene, transcripts in gene_transcripts.items():
            for transcript in transcripts:
                if transcript in csq_transcripts:
                    csq_field = csq_fields[target_field_index]
                    target_transcript_csqs.append( csq_field)
    return target_transcript_csqs


def get_csq_dict(vcf_record, gene_transcripts, csq_field_keys):
    # Needs to handle multi alts
    
    alts = vcf_record.ALT
    csqs = vcf_record.INFO.get("CSQ")
    alleles = alts2csqallele(vcf_record)

    csq_dict = {}

    for csq in csqs:
        transcript_field_name = "RefSeq"
        transcript_field_index = csq_field_keys.index(transcript_field_name)
        
        csq_fields = csq.split("|")
        csq_transcripts = csq_fields[transcript_field_index]
        
        for gene, transcripts in gene_transcripts.items():
            for transcript in transcripts:
                if transcript in csq_transcripts:
                    # it's the right transcript
                    # which alt does it belong to
                    allele_field_name = "Allele"
                    allele_field_index = csq_field_keys.index(allele_field_name)
                    allele = csq_fields[allele_field_index]
                    allele_index = alleles.index(allele)
                    alt = alts[allele_index]
                    csq_dict.setdefault(alt, {})
                    csq_dict[alt][transcript] = csq
    return csq_dict


def generate_consequence_col(csq_field_keys, gene_transcripts):
    def consequence_col(vcf_record):

        csqs = get_csq_field("Consequence", vcf_record, gene_transcripts, csq_field_keys)
        print csqs
        csqs = [csq.split("&") for csq in csqs]
        return csqs

    return consequence_col


def write_row_to_worksheet(worksheet, row):
    '''
    Can't be done with xlwt
    If header is present
        Check against provided fields
    Else
        Write header using provided fields
    
    Find last written line pos and gene
    If same gene
        Write to next line
    else
        Write to next +1 line

    '''
    pass

def alts2csqallele(vcf_record):
    """Convert vcf alts to VEP csq alleles
    
    Args:
        vcf_record (_Record): a pyvcf _Record 
    
    Returns:
        TYPE: An ordered list containing one allele for each alt
    """

    csq_alleles = []
    
    for alt in vcf_record.ALT:
        
        if len(vcf_record.REF) == 1:
            if len(alt) == 1:
            
            # bases in 
            # ref:alt

            # 1:1
                csq_alt = alt
            else:
            # 1:many
                csq_alt = alt[1:]
        else:
            if len(alt) == 1:
            # many:1
                csq_alt == "-"
            else:
            # many:many
                csq_alt == alt[:1]
    
        csq_alleles.append(csq_alt)
    return csq_alleles

def write_row_to_worksheet(ws, row_number, row_data):
    for col_index, value in enumerate(row_data):
        ws.write(row_number, col_index, value)

@begin.start
def main(g_number, vcf_filepath, samplepanels_filepath=None, genes2transcripts_filepath=None, genepanels_filepath=None, panels=None):
    
    # Data files:
    # Which panels are required for each sample
    samplepanels_filepath = samplepanels_filepath or\
                            "/data/gemini/BioinformaticManifest.txt"

    # Which transcript is used for each gene                            
    genes2transcripts_filepath = genes2transcripts_filepath or\
                                 "/data/gemini/genes2transcripts"
    
    # Which genes are in each panel
    genepanels_filepath = genepanels_filepath or\
                                 "/data/gemini/genepanels"
 
    # TO DO - Need to handle reanalysis panels here

    #### Generate sample/panel specific data ####
    panel_IDs = gnum2pids(g_number, samplepanels_filepath)
    panel_names = pids2panelnames(panel_IDs, genepanels_filepath)
    panel_genes = pids2genes(panel_IDs, genepanels_filepath)
    gene_transcripts = genes2transcripts(panel_genes, genes2transcripts_filepath)
    transcript_regions = transcripts2exons(gene_transcripts)
    # Maybe make this into a single dict
  
    #### Generate vcf specific metadata ####
    #csq_fields_index = get_csq_index(vcf_filepath)
    csq_field_keys = get_csq_field_keys(vcf_filepath)

    '''
    # Checking if it looks right - does nothing useful
    for panel in panel_IDs:
        print "\n\nPANEL", panel
        print panel_names[panel]
        for gene in panel_genes[panel]:
            print "\nGene", gene
            for transcript in gene_transcripts[gene]:
                print "Tx", transcript
                for exon_num, region in sorted(transcript_regions[transcript].items()):
                    print "\nExon", exon_num
                    print "%s:%s-%s" % (region["chr"], region["start"],region["end"])
    '''

    # Execute higher order functions to generate col_ functions 
    # These require sample specific data embedded in the function so
    # cannot be generated until this point
    gene_col        = generate_gene_col(gene_transcripts)
    nuc_change_col  = generate_nuc_change_col(csq_field_keys, gene_transcripts)
    aa_change_col   = generate_aa_change_col(csq_field_keys, gene_transcripts)
    consequence_col = generate_consequence_col(csq_field_keys, gene_transcripts)

    # Process vcf records
    with open(vcf_filepath) as vcf_fh:
        vcf_reader = vcf.Reader(vcf_fh)


        # For each column give a column name, and a function which
        # generates the value for the column when given a 
        # pycvf _Record object


        # TODO - Needs a formatting function for each col too
        # alts and transcripts are special cases - their content defines 
        # whether other fields contain multiple values
        column_functions = [("Chromosome", chromosome_col),
                            ("Start", start_col),
                            ("Genomic Ref Allele", ref_col),
                            #("Genomic Alt Allele", alt_col),
                            ("Gene", gene_col),
                            #("Nucleotide pos", nuc_change_col),
                            #("AA Change", aa_change_col),
                            ("Score", qual_score_col),
                            ("AAF", aaf_col),
                            #("Genotype", genotype_col),
                            ("dbsnp", dbsnp_col),
                            #("AF_gemini", af_gemini_col),
                            #("Consequences", consequence_col),
                            ]

        # Set up workbook and sheets

        # For each variant sheet, provide a sheet name and a function 
        # which takes a row and returns the sheet to which it should be 
        # assigned. These will be executed in the order of the list, and
        # the row assigned to the sheet corresponding to the first 
        # function to return True. The final sheet function should
        # always return True to capture any variants not assigned to 
        # other sheets

        # Cell styles
        style_strings = {"bold":"font: bold on",
                         "red":"font: color red",
                        }
        styles = {style_name:xlwt.easyxf(style_strings[style_name]) for style_name in style_strings.keys()}

        # Sheets and filters
        sheet_names = [ ("Test1", True),
                        ("Test2", True),
                        ("Test3", True),
                        ("Test4", True),
                        ("Other", True),
                      ]
        # Make a workbook
        wb = xlwt.Workbook()
        # Make sheets in the workbook
        sheets = {sheet_name[0]:wb.add_sheet(sheet_name[0]) for sheet_name in sheet_names}

        # Need to track the current row in each of the sheets to keep
        # track of position as we write records to different sheets
        current_row_index = {sheet_name[0]:0 for sheet_name in sheet_names}
        
        # For each sheet...
        # (Will need to assign to one sheet for prod)
        for sheet_name, ws in sheets.items():
            
            # ...write the header row
            header_style = styles["red"]
            
            # Each cell is written individually, so go through each col
            # in the row
            for col_index, col in enumerate(column_functions):
                col_name, column_function = col
                ws.write(current_row_index[sheet_name], col_index, col_name, style=header_style)
            current_row_index[sheet_name] += 1
        
        wb.save("test.xls")

        # Now populate the sheets with variants
        for vcf_record in vcf_reader:

            ws = sheets["Test1"]
            
            if record_in_regions(vcf_record, transcript_regions):
                row_data = []
                for c in column_functions:
                    column_header = c[0]
                    column_function = c[1]
                    column_name, column_value = column_header, column_function(vcf_record)
                    row_data.append(column_value)

                for col_index, field in enumerate(row_data):
                    print field

                write_row_to_worksheet(ws, current_row_index[ws.name], row_data)
                current_row_index[ws.name] += 1
                # assign variant to a sheet
                # write it to that sheet
        
        wb.save("test.xls")
#    print csq_fields_index
