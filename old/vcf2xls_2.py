"""Summary
"""
import os
import sys
import MySQLdb
import begin
import vcf
import xlwt
import urllib2



#### Database queries ####

def geminiDB_query(query):
    """Run a query on GeminiDB database
    
    Args:
        query (str): SQL query
    
    Returns:
        list: A list of tuples, each representing one record in the query result
              i.e.[(field1, field2),
                   (field1, field2)...]
    """
    # Database details
    db = MySQLdb.connect(host="mgsrv01",
                         user="easih_ro",
                         db="GeminiDB")

    # Run the query
    result = execute_query(query, db)
    return result


def varfreq_query(query):
    """Run a query on VarFreq database
    
    Args:
        query (str): SQL query
    
    Returns:
        list: A list of tuples, each representing one record in the query result
              i.e.[(field1, field2),
                   (field1, field2)...]
    """
    # Database details
    db = MySQLdb.connect(host="mgsrv01",
                         user="easih_ro",
                         db="VarFreq")

    # Run the query
    result = execute_query(query, db)
    return result


def execute_query(query, db):
    """Run a query on a MySQL database connection
    
    Args:
        query (str): An SQL query
        db (MySQLdb.connect): A connection to a MySQL database
    
    Returns:
        list: A list of tuples, each representing one record in the query result
              i.e.[(field1, field2),
                   (field1, field2)...]
    """
    # Create cursor
    cur = db.cursor()

    # Run query
    cur.execute(query)

    # Get result
    result = []
    for row in cur.fetchall():
        result.append(row)

    db.close()
    return result    


#### Text file data manipulations ####

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
                
                p_name, p_id, p_gene = line_components
                if pid == p_id:
                    panel_name = p_name
                    break

    return panel_name


def panel_name2pid(panel_name, genepanels_filepath):
    if panel_name.startswith("_"):
        pid = panel_name

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


def record_in_regions(vcf_record, regions, flank=0):
    vcf_record_range = record_range(vcf_record)
    
    region_hits = []
    for region, exons in regions.items():
        for exon, coords in sorted(exons.items()):

            # Wrong chrom
            if vcf_record.CHROM != coords["chr"]:
                break
            
            # Right chrom - do coords overlap?
            else:
                region_range = [coords["start"]-flank, coords["end"]+flank]
                
                overlap_length = overlap(region_range, vcf_record_range)
                if overlap_length:
                    print "v", vcf_record_range
                    print "p", region_range
                    print overlap_length
                    region_hits.append(region)
                    
    return region_hits
                   

def record_range(vcf_record):
    alleles = vcf_record.ALT + [vcf_record.REF]
    max_allele_length = max([len(allele) for allele in alleles])
    
    first = vcf_record.POS
    last = vcf_record.POS + max_allele_length

    return [first, last]


def overlap(a, b):
    a = map(int, [x for x in a])
    b = map(int, [x for x in b])
    return max(0, min(a[1], b[1]) - max(a[0], b[0]))


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

# functions should be structured like:
# def field_col(vcf_record, alt, transcript, csq, csq_field_keys)
    # Ignored params
    # del ...
    # do stuff
    # return value_for_col


## vcf field functions
## These only require a vcf_record as input
## Generally DNA level info

def chromosome_col(vcf_record, alt, transcript, csq, csq_field_keys):
    # Ignored params
    del alt, transcript, csq, csq_field_keys
    return vcf_record.CHROM


def start_col(vcf_record, alt, transcript, csq, csq_field_keys):
    # Ignored params
    del alt, transcript, csq, csq_field_keys
    text = vcf_record.POS
    URL = "https://genome.ucsc.edu/cgi-bin/hgTracks?db=hg19&position=chr%s%%3A%d-%d"\
            % (vcf_record.CHROM, vcf_record.POS - 10, vcf_record.POS + 10)
    link = xlwt.Formula('HYPERLINK("%s";"%s")' % (URL, text))

    return link


def ref_col(vcf_record, alt, transcript, csq, csq_field_keys):
    # Ignored params
    del alt, transcript, csq, csq_field_keys
    return vcf_record.REF


def pos_col(vcf_record, alt, transcript, csq, csq_field_keys):
    # Ignored params
    del alt, transcript, csq, csq_field_keys
    return vcf_record.POS


def qual_score_col(vcf_record, alt, transcript, csq, csq_field_keys):
    # Ignored params
    del alt, transcript, csq, csq_field_keys
    return vcf_record.QUAL


def aaf_col(vcf_record, alt, transcript, csq, csq_field_keys):
    # Ignored params
    del alt, transcript, csq, csq_field_keys
    # How to handle multi alts?
    for tag in ["ABHet", "ABHom"]:
        aaf = vcf_record.INFO.get(tag)
        if aaf:
            return aaf


def genotype_col(vcf_record, alt, transcript, csq, csq_field_keys):
    # Ignored params
    del alt, transcript, csq, csq_field_keys
    return [call.data.GT for call in vcf_record.samples]


def dbsnp_col(vcf_record, alt, transcript, csq, csq_field_keys):
    # Ignored params
    del alt, transcript, csq, csq_field_keys
    rsID = vcf_record.ID
    if rsID:
        ID = rsID[2:]  # remove the rs prefix
        URL = "https://www.ncbi.nlm.nih.gov/projects/SNP/snp_ref.cgi?rs=%s" % ID
        rsID = xlwt.Formula('HYPERLINK("%s";"%s")' % (URL, rsID))
    
    return rsID


def af_gemini_col(vcf_record, alt, transcript, csq, csq_field_keys):
    # Ignored params
    del transcript, csq, csq_field_keys
    
    allele_count = gemini_allele_count(vcf_record, alt)

    variant_allele_count = allele_count["hom"]*2 + allele_count["hem"]

    # Ideally would take XX/XY into account but for now assume all diploid
    allele_freq = "%.4f" % (float(variant_allele_count)/(allele_count["total"]*2))
    return allele_freq


def het_hom_gemini_col(vcf_record, alt, transcript, csq, csq_field_keys):
    # Ignored params
    del transcript, csq, csq_field_keys
    
    allele_count = gemini_allele_count(vcf_record, alt)

    return "%d / %d" % (allele_count["hem"], allele_count["hom"])


def af_exac_col(vcf_record, alt, transcript, csq, csq_field_keys):
    ExAC = "ExAC 0.3"

    vids = get_varfreq_vids(vcf_record.CHROM, vcf_record.POS, vcf_record.REF, alt)
    assert len(vids) <= 1, "Multiple vids for one variant in VarFreq"
    
    if vids:
        vid = vids[0]
    
        query = """
                SELECT frequency
                FROM freq
                WHERE vid = '%s'
                AND sid = ( SELECT sid 
                            FROM source 
                            WHERE name = '%s')
                """ % (vid, ExAC)
        
        freq = varfreq_query(query)[0][0]  # Something like [(0.0123,)]
    
    else:
        freq = "0"

    URL = "http://exac.broadinstitute.org/variant/%s-%s-%s-%s" % (vcf_record.CHROM, vcf_record.POS, vcf_record.REF, alt)
    text = freq
    link = xlwt.Formula('HYPERLINK("%s";"%s")' % (URL, text))

    return link


def comment_col(vcf_record, alt, transcript, csq, csq_field_keys):
    pass


def depth_col(vcf_record, alt, transcript, csq, csq_field_keys):
    return [str(call.data.DP) for call in vcf_record.samples]


## csq field functions
## There require a specific vcf_record csq as input
## Generally RNA/protein level info

def gene_col(vcf_record, alt, transcript, csq, csq_field_keys):
    # Ignored params
    del vcf_record, alt, transcript
    field_name = "HGNC"
    gene = get_csq_field(field_name, csq, csq_field_keys)
    return gene


def c_pos_col(vcf_record, alt, transcript, csq, csq_field_keys):
    # Ignored params
    del vcf_record, alt, transcript
    field_name = "HGVSc"
    c_dot = get_csq_field(field_name, csq, csq_field_keys).split(":")[1]
    
    # Find the first letter (excluding c.)
    for index, char in enumerate(c_dot[2:]): # skip the "c." prefix
        if char.isalpha():  
            end_of_pos = index + 2
            break
    
    # Split c_dot into pos and change by splitting at this pos
    c_pos = c_dot[:end_of_pos]
    return c_pos


def c_change_col(vcf_record, alt, transcript, csq, csq_field_keys):
    # Ignored params
    del vcf_record, alt, transcript
    field_name = "HGVSc"
    c_dot = get_csq_field(field_name, csq, csq_field_keys).split(":")[1]
    
    # Find the first letter (excluding c.)
    for index, char in enumerate(c_dot[2:]): # skip the "c." prefix
        if char.isalpha():   
            end_of_pos = index + 2
            break
    
    # Split c_dot into pos and change by splitting at this pos
    c_change = c_dot[end_of_pos:]
    return c_change


def aa_change_col(vcf_record, alt, transcript, csq, csq_field_keys):
    # Ignored params
    del vcf_record, alt, transcript
    field_name = "Amino_acids"
    csq_field = get_csq_field(field_name, csq, csq_field_keys)
    aa_change = csq_field.split("/")
    if len(aa_change) > 1:
        field_name = "Protein_position"
        aa_position = get_csq_field(field_name, csq, csq_field_keys)
        return aa_change[0] + aa_position + aa_change[1]
    else:
        return aa_change


def consequences_col(vcf_record, alt, transcript, csq, csq_field_keys):
    # Ignored params
    del vcf_record, alt, transcript
    field_name = "Consequence"
    consequences = get_csq_field(field_name, csq, csq_field_keys)
    return consequences


def sift_col(vcf_record, alt, transcript, csq, csq_field_keys):
    # Ignored params
    del vcf_record, alt, transcript
    field_name = "SIFT"
    sift = get_csq_field(field_name, csq, csq_field_keys)
    return sift


def polyphen_col(vcf_record, alt, transcript, csq, csq_field_keys):
    # Ignored params
    del vcf_record, alt, transcript
    field_name = "PolyPhen"
    polyphen = get_csq_field(field_name, csq, csq_field_keys)
    return polyphen


def OMIM_col(vcf_record, alt, transcript, csq, csq_field_keys):
    # Use OMIM API to pull inheritance mode info for gene
    
    gene = get_csq_field("HGNC", csq, csq_field_keys)
    print gene
    # Prepare request for gene
    api_request = construct_api_request(gene)

    # Submit request
    api_response = submit_api_request(api_request)
    
    # Extract relevant data from response
    extracted_data = extract_data(api_response, gene)

    # Get inheritance modes
    omim_strings = []
    for IM, disease_list in extracted_data.items():
        IM_string = "%s: %s" % (IM, ", ".join(disease_list))
        omim_strings.append(IM_string)
    omim_string = " | ".join(omim_strings)
    return omim_string


## Other field functions
## These require the specific alt and/or transcript as input

def transcript_col(vcf_record, alt, transcript, csq, csq_field_keys):
    # Ignored params
    del vcf_record
    field_name = "RefSeq"
    csq_transcripts = get_csq_field(field_name, csq, csq_field_keys).split("&")
    assert transcript in csq_transcripts, "Mismatched transcipts!"
    return transcript


def alt_col(vcf_record, alt, transcript, csq, csq_field_keys):
    # Ignored params
    del vcf_record, transcript, csq, csq_field_keys
    return alt.sequence


def af_gnomad_exome_col(vcf_record, alt, transcript, csq, csq_field_keys):
    URL = "http://gnomad.broadinstitute.org/variant/%s-%s-%s-%s" % (vcf_record.CHROM, vcf_record.POS, vcf_record.REF, alt)
    text = "Link"
    link = xlwt.Formula('HYPERLINK("%s";"%s")' % (URL, text))
    return link


def af_gnomad_genome_col(vcf_record, alt, transcript, csq, csq_field_keys):
    URL = "http://gnomad.broadinstitute.org/variant/%s-%s-%s-%s" % (vcf_record.CHROM, vcf_record.POS, vcf_record.REF, alt)
    text = "Link"
    link = xlwt.Formula('HYPERLINK("%s";"%s")' % (URL, text))
    return link


########################

# Non-variant sheets
def add_summary_sheet(wb, styles):
    
    ws = wb.add_sheet("Summary")

    # Individual cells

                # Value             # Style              # Row Col
    cells = [   ("Gemini ID:",      styles["bold"],         0,  0),
                #("Value1",          styles["red_text"],     0,  1),
                
                ("GM number:",      styles["bold"],         1,  0),
                #("Value2",          styles["red_text"],     1,  1),
                
                ("Name:",           styles["bold"],         2,  0),
                #("Value3",          styles["red_text"],     2,  1),
            
                ("Inferred Gender:",styles["bold"],         0,  2),
                #("Value4",          styles["red_text"],     0,  3),
                
                ("SRY present:",    styles["bold"],         1,  2),
                #("Value5",          styles["red_text"],     1,  3),

                ("Panel Coverage:", styles["bold"],         0,  4),
                #("Value6",          styles["red_text"],     0,  5),

                ("Panels:",         styles["bold"],         0,  6),
                #("Value7",          styles["red_text"],     0,  7),

                ("Report Text:",    styles["bold"],         6,  0),
                #("Value8",          styles["red_text"],     6,  1),

                ("Bics/Seq QC",     styles["bold"],         38,  0),
                ("Reads",           styles["bold"],         39,  0),
                #("Value9",          styles["red_text"],     39,  1),
                ("Usable Reads",    styles["bold"],         40,  0),
                #("Value10",         styles["red_text"],     40,  1),

            ]


    for cell in cells:
        value, style, row, col = cell
        ws.write(row, col, value, style)

    # Tables

    # Phenotype table
    ws = add_table(ws, (8,1), (5,4), ["Phenotype"], [(0,3)], styles)

    # panels table
    column_headers = ["Panels", "Excel File", "Comments", "Analysis By", "Date", "Checked By", "Date"]
    ws = add_table(ws, (15,1), (3,8), column_headers, [(2,3)], styles)
    
    # Sanger table extra header
    ws = add_table(ws, (20,1), (1,7), ["Sanger sequencing confirmation"], [(0,6)], styles)

    # Sanger table
    column_headers = ["Gene","NM_#","coordinate","cDNA","protein change","WS#","confirmed(y/n)"]
    ws = add_table(ws, (21,1), (4,7), column_headers, [], styles)

    # GEM comments table
    column_headers = ["GEM comments summary", "date"]
    ws = add_table(ws, (27,1), (5,5), column_headers, [(0,1),(2,4)], styles)
    
    return wb


def add_geneqc_sheet(wb, sample, styles, gene_transcripts):
    
    ws = wb.add_sheet("Gene_QC")

    columns = [ "Name","Transcript","Region length","Min depth",
                "Missing","1-5x","6-9x","10-19x","20+x","20+x %"]
    
    row_index = 0

    # Write headers
    for col_index, column_header in enumerate(columns):
        ws.write(row_index, col_index, column_header, style=styles["bold"])
    row_index += 1

    for gene, transcripts in sorted(gene_transcripts.items()):
        for transcript in transcripts:
            gene_performance_data = get_gene_performance(sample, transcript)
            
            # Set style rules here
            # Default
            row_style = styles["plain_text"]
            # Exceptions
            if float(gene_performance_data["Gene"]["20+x %"]) < 100.0:
                row_style = styles["red_text"]

            # Write data rows
            for col_index, column_header in enumerate(columns):
                value = gene_performance_data["Gene"][column_header]
                ws.write(row_index, col_index, value, row_style)
            row_index += 1
        row_index += 1

    ws = freeze_header(ws)
    return wb


def add_exonqc_sheet(wb, sample, styles, gene_transcripts):
    
    ws = wb.add_sheet("Exon_QC")

    columns = [ "Name","Transcript","Position","Min depth","Max depth",
                "Mean depth", "Exp mean depth", "Missing","1-5x","6-9x",
                "10-19x", "CNV type", "CNV region", "CNV BF", "CNV Confidence"]
    
    row_index = 0

    # Write headers
    for col_index, column_header in enumerate(columns):
        ws.write(row_index, col_index, column_header, style=styles["bold"])
    row_index += 1

    # Get data
    for gene, transcripts in sorted(gene_transcripts.items()):
        for transcript in transcripts:
            gene_performance_data = get_gene_performance(sample, transcript)
            
            # Write data rows
            for exon in sorted(gene_performance_data["Exon"]):
                
                # Set style rules here
                # Default
                row_style = styles["plain_text"]
                # Exceptions
                if int(gene_performance_data["Exon"][exon]["Min depth"]) < 20:
                    row_style = styles["red_text"]

                for col_index, column_header in enumerate(columns):
                    value = gene_performance_data["Exon"][exon].get(column_header, "")
                    ws.write(row_index, col_index, value, row_style)
                row_index += 1
            row_index += 1

    ws = freeze_header(ws)

    return wb


def add_table(ws, position, dimensions, col_headers, merged_cols, styles):
    """Add a table to a sheet
    
    Args:
        ws (xlwt.Workbook.Worksheet): an xl worksheet
        position (tuple): (row, col) of the upper left corner of the table (0 indexed)
        dimensions (tuple): (num_rows, num_cols) in the table (prior to column merges)
        col_headers (tuple): The column headers (as strings) for the table
        merged_cols (tuple): First and last column in each set of cols to be merged (i.e. ((1,3),(6,7)) will merge cols 1-3 into one col, and cols 6-7 into one col)
        styles (dict): Styles to be applied to cells
    
    Returns:
        xlwt.Workbook.Worksheet: The worksheet containing the table
    """

    header_style = styles["summary_table_header"]
    cell_style = styles["table_cell"]
    table_edge = styles["table_merge_fix"]  # See last loop in func

    first_row, first_col = position
    num_rows, num_cols = dimensions
    rows = [first_row + n for n in range(0, num_rows)]
    cols = [first_col + n for n in range(0, num_cols)]

    # No point writing values or formatting cells which are going to merge
    # with others, as this info is lost for all but top left cell of merge
    cols_lost_to_merge = []
    for merge_range in merged_cols:
        for col in range(first_col+merge_range[0]+1, first_col+merge_range[1]+1):
            cols_lost_to_merge.append(col)

    # Check that we have enough col headers
    num_headers_required = len(cols) - len(cols_lost_to_merge)
    assert len(col_headers) == num_headers_required, "%d, %d" % (len(col_headers), num_headers_required)

    
    header_row = rows[0]
    col_index = -1  # Tracking cols excluding those which will merge

    # Merge cells
    for merge_range in merged_cols:
        left_col, right_col = [x + first_col for x in merge_range]
        for row in rows:
            ws.merge(row,row,left_col,right_col)

    # Write values and formats
    for col in cols:
        if col not in cols_lost_to_merge:
            col_index += 1
            col_header = col_headers[col_index]

            # Write header
            ws.write(header_row, col, col_header, header_style)

            for row in rows[1:]:
                ws.write(row, col, style=cell_style)

    # Merging cells causes the right border on the last col of the table
    # to be lost. Seems to be a xlwt bug. 
    # Workaround:

    next_col = cols[-1]+1
    for row in rows:
        ws.write(row, next_col, style=table_edge)

    return ws


def freeze_header(ws):
    ws.set_panes_frozen(True)
    ws.set_horz_split_pos(1)
    return ws


def get_varfreq_vids(chrom, pos, ref, alt):
    query = """
            SELECT vid
            FROM variant
            WHERE chr = '%s'
            AND pos = '%s'
            AND ref = '%s'
            AND alt = '%s'
            """ % (chrom, pos, ref, alt)

    vids = [x[0] for x in varfreq_query(query)]
    return vids


def gemini_allele_count(vcf_record, alt):
    gemini_prid = "1"

    # Get variant ID(s)
    query = """ SELECT vid
                FROM variant
                WHERE chr = '%s'
                AND pos = '%s'
                AND ref = '%s'
                AND alt = '%s'
            """ % (vcf_record.CHROM, vcf_record.POS, vcf_record.REF, alt)
    
    vids = [line[0] for line in geminiDB_query(query)]
    
    # Some variants have multiple vids for some reason
    # The same sample may have the same variant against both vids
    # Therefore we generate a variant count for each vid, and take the 
    # highest

    variants = []
    
    for vid in vids:
    
        query = """ SELECT s_v.sid, s_v.allele_count
                    FROM sample_variant as s_v
                    INNER JOIN sample as sam
                    ON s_v.sid = sam.sid
                    WHERE s_v.vid = '%s'
                    AND sam.prid = '%s'
                """ % (vid, gemini_prid)
    
        result = geminiDB_query(query)

        # Keep the longest list
        if len(result) > len(variants):
            variants = result
    
    # Count hem and hom samples
    genotype_count = {  "hem":   0,
                        "hom":   0,
                     }

    for variant in variants:
        allele_count = variant[1]
        if allele_count == 1:
            genotype_count["hem"] += 1
        elif allele_count == 2:
            genotype_count["hom"] += 1
        else:
            assert False, "Allele count is %s!" % allele_count
    
    # Also need to know how many samples total to calc ratios
    
    query = """ SELECT count(sid)
                FROM sample
                WHERE prid = '%s'
            """ % (gemini_prid)
    
    gemini_sample_count = geminiDB_query(query)[0][0]
    genotype_count["total"] = gemini_sample_count

    return genotype_count


# CSQ stuff

def get_csq_field(field_name, csq, csq_field_keys):

    target_field_index = csq_field_keys.index(field_name)
    csq_fields = csq.split("|")
    csq_target_value = csq_fields[target_field_index]    
    return csq_target_value


def get_csq_dict(vcf_record, gene_transcripts, csq_field_keys):
    
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


def alts2csqallele(vcf_record):
    """Convert vcf alts to VEP csq alleles
    
    Args:
        vcf_record (_Record): a pyvcf _Record 
    
    Returns:
        TYPE: An ordered list containing one allele for each alt
    """

    csq_alleles = []
    
    for alt in vcf_record.ALT:
        
        alt_seq = alt.sequence

        if len(vcf_record.REF) == 1:
            if len(alt_seq) == 1:
            
            # bases in 
            # ref:alt

            # 1:1
                csq_alt = alt_seq
            else:
            # 1:many
                csq_alt = alt_seq[1:]
        else:
            if len(alt_seq) == 1:
            # many:1
                csq_alt = "-"
            else:
            # many:many
                csq_alt = alt_seq[:1]
    
        csq_alleles.append(csq_alt)
    return csq_alleles


def write_row_to_worksheet(ws, row_number, row_data, style):
    for col_index, value in enumerate(row_data):
        # print "Row:\t%d\tCol:\t%d\tValue:\t%s" % (row_number, col_index, value)
        ws.write(row_number, col_index, value, style)



# Decon Stuff
   
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

    return decon_data


def get_decon_results_filepath(sample):
    # TO DO
    return "/data/projects/matt/cnv/DECoN_results/170405_170407_NS500192/170405_170407_NS500192_all.txt"


def update_exon_with_decon(exon_data, decon_data_dict):
    # Set default values to be used if no CNV found
    exon_data["CNV region"] = ""
    exon_data["CNV type"] = ""
    exon_data["CNV BF"] = ""
    exon_data["CNV Confidence"] = ""

    # Search for a CNV in the exon
    for cnv_id, cnv_data in decon_data_dict.items():
        # Same sample and chrom?
        if (cnv_data["Sample"] == exon_data["Sample"]) and (cnv_data["Chromosome"] == exon_data["chr"]):
            # Do coords of exon overlap those of cnv?
            region_overlap = overlap((exon_data["start"], exon_data["end"]), (cnv_data["Start"], cnv_data["End"]))
            
            if region_overlap:
                cnv_region = "%s:%s-%s" %  (cnv_data["Chromosome"], cnv_data["Start"], cnv_data["End"])
                
                exon_data["CNV type"] = cnv_data["CNV.type"]
                exon_data["CNV region"] = cnv_region
                exon_data["CNV BF"] = cnv_data["BF"]
                exon_data["CNV Confidence"] = decon_bayes_interpretation(cnv_data["BF"])
                break

    return exon_data


def decon_bayes_interpretation(bayes_factor):
    # Based on Kass and Raftery (1995)

    bayes_factor = float(bayes_factor)

    if bayes_factor < 1:
        return "Insignificant"
    elif 1 <= bayes_factor <= 3:
        return "Weak"
    elif 3 < bayes_factor <= 20:
        return "Positive"
    elif 20 < bayes_factor <= 150:
        return "Strong"
    elif bayes_factor > 150:
        return "Very Strong"
    else:
        return "ERROR"


# Sequencing QC stuff

def get_gene_performance(sample, transcript):
    
    # Query the gemini database for performance data
    query = """ SELECT  sam.name as sample, 
                        ge.name as HGNC, 
                        r2g.exon_nr as exon, 
                        ge.refseq as RefSeq, 
                        reg.*, 
                        cov.min, 
                        cov.mean, 
                        cov.max, 
                        cov.missing, 
                        cov.1to5, 
                        cov.6to9, 
                        cov.10to19 
                FROM    coverage as cov 
                        inner join sample as sam on sam.sid = cov.sid 
                        inner join region2gene as r2g on cov.rid = r2g.rid 
                        inner join gene as ge on ge.gid = r2g.gid 
                        inner join region as reg on reg.rid = r2g.rid 
                WHERE   ge.refseq like '%%%s%%' and sam.name = '%s' 
            """ % (transcript, sample)
    result = geminiDB_query(query)

    # Put output into a dict
    performance_data = {"Gene":{},  # Gene level info in this dict
                        "Exon":{},  # Exon level info in this dict
                        "Sample":None,  # Just the sample ID?
                       }

    decon_filepath = get_decon_results_filepath(sample)
    decon_data_dict = decon_dict(decon_filepath)

    first_record = True
    for record in result:

        # We get one record per exon
        sample, HGNC, exon, refseq, rid, chr_, start, end, min_, mean,\
        max_, missing, oneto5, sixto9, tento19 = record

        if first_record:
            performance_data["Sample"] = sample
            performance_data["Gene"]["Name"] = HGNC
            performance_data["Gene"]["Transcript"] = transcript
            first_record = False

        # Pull values straight from query
        performance_data["Exon"][exon] = { "Sample":sample,
                                           "Gene":HGNC,
                                           "Transcript":transcript,
                                           "Exon":exon,
                                           "rid":rid,
                                           "chr":chr_,
                                           "start":start,
                                           "end":end,
                                           "Min depth":min_,
                                           "Mean depth":mean,
                                           "Max depth":max_,
                                           "0x":missing,
                                           "1-5x":oneto5,
                                           "6-9x":sixto9,
                                           "10-19x":tento19,
                                           }

        # Calc values not in query
        performance_data["Exon"][exon]["Name"] = calc_exon_name(performance_data["Exon"][exon])
        performance_data["Exon"][exon]["Position"] = calc_exon_position(performance_data["Exon"][exon])
        performance_data["Exon"][exon]["Exp mean depth"] = calc_exon_exp_mean_depth(performance_data["Exon"][exon])

        # Get DECoN dosage calls for exons - adds multiple fields
        performance_data["Exon"][exon] = update_exon_with_decon(performance_data["Exon"][exon], decon_data_dict)
        

    # Once the full exon dict is built we then generate values for each
    # gene(# missing bases etc) by using the dict

    # The order of these calls matters - some are dependent on the
    # output of preceeding calls being present in the performance data dict
    performance_data["Gene"]["Region length"] = calc_region_length(performance_data)
    performance_data["Gene"]["Min depth"] = calc_gene_min_depth(performance_data)
    performance_data["Gene"]["Missing"] = calc_bases_in_depth_bin(performance_data, "0x")
    performance_data["Gene"]["1-5x"] = calc_bases_in_depth_bin(performance_data, "1-5x")
    performance_data["Gene"]["6-9x"] = calc_bases_in_depth_bin(performance_data, "6-9x")
    performance_data["Gene"]["10-19x"] = calc_bases_in_depth_bin(performance_data, "10-19x")
    performance_data["Gene"]["20+x"] = calc_gene_20x_len(performance_data)
    performance_data["Gene"]["20+x %"] = calc_gene_20x_perc(performance_data)    
    
    return performance_data


def get_gene_performance2(sample, transcript):
    
    # Query the gemini database for performance data
    query = """ SELECT  sam.name as sample, 
                        ge.name as HGNC, 
                        r2g.exon_nr as exon, 
                        ge.refseq as RefSeq, 
                        reg.*, 
                        cov.min, 
                        cov.mean, 
                        cov.max, 
                        cov.missing, 
                        cov.1to5, 
                        cov.6to9, 
                        cov.10to19 
                FROM    coverage as cov 
                        inner join sample as sam on sam.sid = cov.sid 
                        inner join region2gene as r2g on cov.rid = r2g.rid 
                        inner join gene as ge on ge.gid = r2g.gid 
                        inner join region as reg on reg.rid = r2g.rid 
                WHERE   ge.refseq like '%%%s%%' and sam.name = '%s' 
            """ % (transcript, sample)
    result = geminiDB_query(query)

    # Put output into a dict
    performance_data = {"Gene":{},  # Gene level info in this dict
                        "Exon":{},  # Exon level info in this dict
                        "Sample":None,  # Just the sample ID?
                       }

    decon_filepath = get_decon_results_filepath(sample)
    decon_data_dict = decon_dict(decon_filepath)

    first_record = True
    for record in result:

        # We get one record per exon
        sample, HGNC, exon, refseq, rid, chr_, start, end, min_, mean,\
        max_, missing, oneto5, sixto9, tento19 = record

        if first_record:
            performance_data["Sample"] = sample
            performance_data["Gene"]["Name"] = HGNC
            performance_data["Gene"]["Transcript"] = transcript
            first_record = False

        # Pull values straight from query
        performance_data["Exon"][exon] = { "Sample":sample,
                                           "Gene":HGNC,
                                           "Transcript":transcript,
                                           "Exon":exon,
                                           "rid":rid,
                                           "chr":chr_,
                                           "start":start,
                                           "end":end,
                                           "Min depth":min_,
                                           "Mean depth":mean,
                                           "Max depth":max_,
                                           "0x":missing,
                                           "1-5x":oneto5,
                                           "6-9x":sixto9,
                                           "10-19x":tento19,
                                           }

        # Calc values not in query
        performance_data["Exon"][exon]["Name"] = calc_exon_name(performance_data["Exon"][exon])
        performance_data["Exon"][exon]["Position"] = calc_exon_position(performance_data["Exon"][exon])
        performance_data["Exon"][exon]["Exp mean depth"] = calc_exon_exp_mean_depth(performance_data["Exon"][exon])

        # Get DECoN dosage calls for exons - adds multiple fields
        performance_data["Exon"][exon] = update_exon_with_decon(performance_data["Exon"][exon], decon_data_dict)
        

    # Once the full exon dict is built we then generate values for each
    # gene(# missing bases etc) by using the dict

    # The order of these calls matters - some are dependent on the
    # output of preceeding calls being present in the performance data dict
    performance_data["Gene"]["Region length"] = calc_region_length(performance_data)
    performance_data["Gene"]["Min depth"] = calc_gene_min_depth(performance_data)
    performance_data["Gene"]["Missing"] = calc_bases_in_depth_bin(performance_data, "0x")
    performance_data["Gene"]["1-5x"] = calc_bases_in_depth_bin(performance_data, "1-5x")
    performance_data["Gene"]["6-9x"] = calc_bases_in_depth_bin(performance_data, "6-9x")
    performance_data["Gene"]["10-19x"] = calc_bases_in_depth_bin(performance_data, "10-19x")
    performance_data["Gene"]["20+x"] = calc_gene_20x_len(performance_data)
    performance_data["Gene"]["20+x %"] = calc_gene_20x_perc(performance_data)    
    
    return performance_data


def calc_region_length(performance_data):
    region_length = sum([exon["end"] - exon["start"] + 1for exon in performance_data["Exon"].values()])
    return region_length


def calc_gene_min_depth(performance_data):
    gene_min_depth = float("inf")

    for exon in performance_data["Exon"].values():
        exon_min_depth = exon["Min depth"]
        if (gene_min_depth == None) or (gene_min_depth > exon_min_depth):
            gene_min_depth = exon_min_depth

    if gene_min_depth == float("inf"):
        gene_min_depth = None 

    return gene_min_depth


def calc_bases_in_depth_bin(performance_data, depth_bin):
    
    bases_in_depth_bin = 0

    for exon in performance_data["Exon"].values():
        regions = exon[depth_bin].split(",")
        for region in regions:
            if region:
            # e.g. 7:117232563-117232564
                start, end = map(int, region.split(":")[1].split("-"))
                bases_in_depth_bin += end-start+1

    return bases_in_depth_bin


def calc_gene_20x_len(performance_data):
    # No 20x+ bin in database, so need to calc this using other bins

    total_len = performance_data["Gene"]["Region length"]
    under_20x_len = 0

    for exon in performance_data["Exon"].values():
        for depth in ["0x", "1-5x", "6-9x", "10-19x"]:
            low_depth_regions = exon[depth].split(",")
            for region in low_depth_regions:
                if region:
                # e.g. 7:117232563-117232564
                    start, end = map(int, region.split(":")[1].split("-"))
                    under_20x_len += end-start+1
    len_20x = total_len - under_20x_len

    return len_20x


def calc_gene_20x_perc(performance_data):
   len_total = performance_data["Gene"]["Region length"]
   len_20x = performance_data["Gene"]["20+x"]
   perc_20x = float(len_20x)/len_total *100
   return round(perc_20x, 2)


def calc_exon_name(exon_data):
    exon_number = exon_data["Exon"]
    gene_name = exon_data["Gene"]
    exon_name = "_".join(map(str ,[gene_name, exon_number]))
    return exon_name


def calc_exon_position(exon_data):
    exon_chrom = exon_data["chr"]
    exon_start = exon_data["start"]
    exon_end = exon_data["end"]
    position = "%s:%d-%d" % (exon_chrom, exon_start, exon_end)
    return position


def calc_exon_exp_mean_depth(exon):
    return "TODO"


# Variant sheet filters

def generate_variant_type_filter(variant_types):
    def variant_type_filter(csq, csq_field_keys):
        field_name = "Consequence"
        consequences = get_csq_field(field_name, csq, csq_field_keys)
        print consequences
        for variant_type in variant_types:
            print "Searching for:\n%s\nin\n%s\n" % (variant_type, consequences)
            if variant_type in consequences:
                print True
                return True
        else:
            print False
            return False
    return variant_type_filter

stop_gained_filter = generate_variant_type_filter(["stop_gained_variant"])
frameshift_filter =  generate_variant_type_filter(["frameshift_variant"])
splice_filter =  generate_variant_type_filter(["splice_region_variant"])
missense_filter = generate_variant_type_filter(["missense_variant"])
synonymous_filter = generate_variant_type_filter(["synonymous_variant"])
# Other_filter should always contain "" to act as a catch-all for variants
# which pass all other filters
other_filter = generate_variant_type_filter([""])


#### OMIM Stuff ####

def construct_api_request(search_term):
    '''
    Assemble an OMIM API URL for a specific search term
    '''
    omim_api_url = "http://api.omim.org/api"
    response_format = "python"
    handler = "entry"
    action = "search"
    API_key = "apiKey=YwfXGOC4Sim_82R9qt8Xzw"
    limit_results = 5
    search_params = "search=" + str(search_term)
    include = "all"

    api_request = omim_api_url + "/" + handler + "/" + action + "?" + search_params  + "&limit=" + str(limit_results) +  "&format=" + response_format + "&include=" + include + "&" + API_key
    return api_request


def submit_api_request(api_request):
    '''
    Submit an api request to the OMIM API and return the response
    '''
    print api_request
    api_response = urllib2.urlopen(api_request)
    response_data = eval(api_response.read())
    return response_data


def extract_data(data, query_HUGO):
    '''
    Extract entries for a specific HGNC gene name from returned OMIM API data
    '''

    output = []
    entries = data["omim"]["searchResponse"]["entryList"]

    for gene_entry in entries:

        # Only process entries for the correct gene
        result_HUGO = gene_entry["entry"]["titles"]["preferredTitle"].split("; ")[-1]
        if result_HUGO == query_HUGO:

            # Gene level data - one per entry
            gene_mimNumber = gene_entry["entry"]["mimNumber"]
            gene_name = gene_entry["entry"]["titles"]["preferredTitle"]
            gene_symbol = gene_entry["entry"]["titles"]["preferredTitle"].split("; ")[-1]

            # Phenotype level data - many per gene/entry
            # .get used to return the value if it exists, otherwise None
            for phenotype in gene_entry['entry']["geneMap"]["phenotypeMapList"]:
                phenotype_name = phenotype.get("phenotypeMap", {}).get("phenotype")
                phenotype_mimNumber = phenotype.get("phenotypeMap", {}).get("phenotypeMimNumber")
                phenotype_inheritance_mode = phenotype.get("phenotypeMap", {}).get("phenotypeInheritance")

                # Additional API call is needed to get phenotype symbol...
                '''
                if phenotype_mimNumber:
                    phenotype_api_request = construct_api_request(phenotype_mimNumber)
                    phenotype_data = submit_api_request(phenotype_api_request)
                    entries = phenotype_data["omim"]["searchResponse"]["entryList"]
                    # Searches for mims always return the mim as result #1 so use entries[0]
                    phenotype_symbol = entries[0]["entry"]["titles"]["preferredTitle"].split("; ")[-1]
                else:
                    phenotype_symbol = None
                '''

                # Don't return 'nondiseases' ' [ ' 
                # or susceptibility for multifactorial diseases ' { '
                if phenotype_name[0] not in ["[", "{"]:
                    output.append([phenotype_name,
                                   #phenotype_symbol,
                                   #phenotype_mimNumber,
                                   phenotype_inheritance_mode,
                                   #gene_name,
                                   #gene_mimNumber,
                                   ])

    inheritance_dict = {}
    for disease, mode in output:
        inheritance_dict.setdefault(mode, []).append(disease)
    return inheritance_dict


@begin.start
def main(g_number, vcf_filepath, samplepanels_filepath=None, genes2transcripts_filepath=None, genepanels_filepath=None, panels=None):
    
    # Data files to tell us...
    # Which panels are required for each sample
    samplepanels_filepath = samplepanels_filepath or\
                            "/data/gemini/BioinformaticManifest.txt"

    # Which transcript is used for each gene                            
    genes2transcripts_filepath = genes2transcripts_filepath or\
                                 "/data/gemini/genes2transcripts"
    
    # Which genes are in each panel
    genepanels_filepath = genepanels_filepath or\
                                 "/data/gemini/genepanels"

    # Variants within exons +/- exon_flank will be included in report
    exon_flank = 30
 
    #### Generate sample/panel specific data ####
    if panels:
        panel_IDs = panel_names2pids(panels)
    else:
        panel_IDs = gnum2pids(g_number, samplepanels_filepath)
    panel_names = pids2panelnames(panel_IDs, genepanels_filepath)
    panel_genes = pids2genes(panel_IDs, genepanels_filepath)
    gene_transcripts = genes2transcripts(panel_genes, genes2transcripts_filepath)
    transcript_regions = transcripts2exons(gene_transcripts)
  
    #### Generate vcf specific metadata ####
    csq_field_keys = get_csq_field_keys(vcf_filepath)

    # Here we specify column headers, functions to generate values,
    # and col widths (in chars)
    # Cols will appear in the same order in the xls report
    column_functions =     [("Gene", gene_col, 8),
                            ("Transcript", transcript_col, 16),
                            ("Chr", chromosome_col, 4),
                            ("Start", start_col, 11),
                            ("Ref", ref_col, 5),
                            ("Alt", alt_col, 5),
                            ("Nucleotide pos", c_pos_col, 12),
                            ("Nucleotide change", c_change_col, 12),
                            ("AA Change", aa_change_col, 10),
                            ("Score", qual_score_col, 9),
                            ("Depth", depth_col, 7),
                            ("AAF", aaf_col, 6),
                            ("Genotype", genotype_col, 9),
                            ("dbsnp", dbsnp_col, 12),
                            ("Polyphen", polyphen_col, 14),
                            ("SIFT", sift_col, 14),
                            ("AF\nGemini", af_gemini_col, 9),
                            ("het/hom\nGemini", het_hom_gemini_col, 12),
                            ("AF\nExAC", af_exac_col, 12),
                            ("AF\nGnomad(E)", af_gnomad_exome_col, 12),
                            ("AF\nGnomad(G)", af_gnomad_genome_col, 12),
                            ("Consequences", consequences_col, 30),
                            ("OMIM", OMIM_col, 30),
                            ]

    column_names = [x[0] for x in column_functions]
    column_widths = [x[2] for x in column_functions]
    

    # Set up the xls workbook
    with open(vcf_filepath) as vcf_fh:
        vcf_reader = vcf.Reader(vcf_fh)
 
 
        # Excel cell formatting - format strings are like css
        style_strings = {"header":"font: bold on;\
                                   align: wrap yes;",
                         "bold":"font: bold on",
                         "plain_text":"font: color black",
                         "red_text":"font: color red",
                         "summary_table_header": "pattern: pattern solid, fore_colour cyan_ega;\
                                                  font: color black;\
                                                  borders: left thin, right thin, top thin, bottom thin;",
                         "table_cell":           "borders: left thin, right thin, top thin, bottom thin;",
                         "table_merge_fix":      "borders: left thin;",
                         "hyperlink": "font: color blue",
                        }

        # Use the strings to generate easyxf objects
        styles = {style_name:xlwt.easyxf(style_strings[style_name]) \
                    for style_name in style_strings.keys()}

        
        # Sheets and filters
        # For each variant sheet, provide a sheet name and a function 
        # which takes a row and returns true if the row should be 
        # assigned to the sheet. 
        # These will be executed in the order of the list, and
        # the row is assigned to the sheet corresponding to the first 
        # function to return True. The final sheet function should
        # always return True to capture any variants not assigned to 
        # other sheets

        sheets = [      ("stop_gained", stop_gained_filter),
                        ("frameshift_variant", frameshift_filter),
                        ("splice", splice_filter),
                        ("missense", missense_filter),
                        ("other", other_filter),
                 ]

        # Make a workbook
        wb = xlwt.Workbook()

        # Add the non-variant sheets
        wb = add_summary_sheet(wb, styles)  
        wb = add_geneqc_sheet(wb, g_number, styles, gene_transcripts)
        wb = add_exonqc_sheet(wb, g_number, styles, gene_transcripts)

        # Add variant sheets
        variant_sheets = {sheet[0]:wb.add_sheet(sheet[0]) for sheet in sheets}

        # Track the current row in each of the sheets to enable 
        # us to write to the correct line as we write records to 
        # different sheets
        current_row_index = {sheet[0]:0 for sheet in sheets}
        
        # Add variant sheet headers
        header_style = styles["header"]
        for sheet_name, ws in variant_sheets.items():
                        
            # Each cell is written individually, so go through each col
            # in the row
            for col_index, col in enumerate(column_functions):
                col_name, column_function, _ = col
                row_index = current_row_index[sheet_name]
                ws.write(row_index, col_index, col_name, style=header_style)
            freeze_header(ws)  # Prevent the header row from scrolling off screen
            current_row_index[sheet_name] += 1
        
            # Set col widths - 1 char ~= 256 units
            for index, col_width in enumerate(column_widths):
                ws.col(index).width = 256 * col_width

            # Set header row height - doesn't affect individial cells in the row
            ws.row(0).height_mismatch = True  # Unlocks height from cell font height
            ws.row(0).height = 256 * 3

        wb.save("test.xls")


        # Now populate the sheets with variants
        for vcf_record in vcf_reader:

            # Get a list of alt alleles in the record
            alts = vcf_record.ALT

            # Get IDs of all transcripts of interest which the variant hits
            transcripts = record_in_regions(vcf_record, transcript_regions, exon_flank)

            # Split csqs according to alt and transcript
            csq_dict = get_csq_dict(vcf_record, gene_transcripts, csq_field_keys)

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
                    if float(row_data[column_names.index("AF\nGemini")]) < 0.01:
                        row_style = styles["red_text"]

                    # Assign variant to a sheet using filter functions
                    for sheet in sheets:
                        sheet_name, sheet_filter = sheet
                        if sheet_filter(csq, csq_field_keys):
                            ws = variant_sheets[sheet_name]
                            write_row_to_worksheet(ws, current_row_index[ws.name], row_data, row_style)
                            current_row_index[ws.name] += 1
                            break

                    

    wb.save("test.xls")
