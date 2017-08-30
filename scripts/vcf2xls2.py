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

def get_output_filepath(vcf_filepath):
    
    # If a filename is not specified we will generate an output filepath
    # in the same dir as the input vcf, with a filename like:
    # G001234_v3.xls

    sample = os.path.basename(vcf_filepath).split(".")[0]
    directory = "/data/projects/matt/vcf2xls_py/test_data/" # debug os.path.dirname(vcf_filepath)
    xls_filename = "%s_v3.xls" % sample
    output_filepath = os.path.join(directory, xls_filename)

    # If the default filename already exists, then add a numerical
    # suffix, increment by 1 until a unique name is found, like:
    # G001234_v3_1.xls  
    i = 0
    while os.path.exists(output_filepath):
        i += 1
        suffix = "_%d" % i
        xls_filename = "%s_v3%s.xls" % (sample, suffix)
        output_filepath = os.path.join(directory, xls_filename)

    return os.path.abspath(output_filepath)

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

def genes2regions(genes, target_bed):
    # Gene -> Tx -> Exon -> Region
    gene_regions = {}
    with open(target_bed) as bed_fh:
        for line in bed_fh:
            chrom, start, end, region_id = line.strip().split("\t")
            gene, enst, refseq_prefix, refseq_suffix, exon = region_id.split("_")
            refseq = refseq_prefix+refseq_suffix
            if gene in genes:
                gene_regions.setdefault(gene, {})
                gene_regions[gene].setdefault(refseq, {})
                gene_regions[gene][refseq][int(exon)] = {"chrom":chrom,
                                                         "start":int(start),
                                                         "end":int(end)}
    return gene_regions

# Database interactions

def execute_query(query, db):
    """Execute a query on a MySQL database connection
    
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
    

def geminiDB_query(query):
    """Query the GeminiDB database
    
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


def gemini_allele_count(chrom, pos, ref, alt):
    """Return the Gemini allele count for a specified variant
    
    Args:
        chrom (str): Chromosome name (without chr prefix)
        pos (str): Position of the variant within the chromosome
        ref (str): The reference allele of the variant
        alt (str): The alternate allele of the variant
    
    Returns:
        dict: Counts of observed genotypes:
                genotype_count{"hem": int
                               "hom": int
                               "total": int
                              }
    """
    
    # GeminiDB database project id for Gemini samples
    gemini_prid = "1"

    # Get variant ID(s)
    query = """ SELECT vid
                FROM variant
                WHERE chr = '%s'
                AND pos = '%s'
                AND ref = '%s'
                AND alt = '%s'
            """ % (chrom, pos, ref, alt)
    
    vids = [line[0] for line in geminiDB_query(query)]
    
    # Some variants have multiple vids for some reason
    # The same sample may have the same variant against both vids
    # Therefore we generate a variant count for each vid, and take the 
    # highest

    variants = []
    
    for vid in vids:
        # Generate a variant count
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


def af_gemini_col(variantrecord, alt, csq, csq_field_keys):
    
    allele_count = gemini_allele_count(variantrecord.chrom, variantrecord.pos, variantrecord.ref, alt)

    variant_allele_count = allele_count["hom"]*2 + allele_count["hem"]

    # Ideally would take XX/XY into account but for now assume all diploid
    allele_freq = "%.4f" % (float(variant_allele_count)/(allele_count["total"]*2))
    return allele_freq


def het_hom_gemini_col(variantrecord, alt, csq, csq_field_keys):
    
    allele_count = gemini_allele_count(variantrecord.chrom, variantrecord.pos, variantrecord.ref, alt)

    return "%d / %d" % (allele_count["hem"], allele_count["hom"])

def varfreq_query(query):
    """Query the VarFreq database
    
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

def get_varfreq_vids(chrom, pos, ref, alt):
    """Get VarFreq database variant ids for a specified variant
    
    Args:
        chrom (str): Chromosome name (without chr prefix)
        pos (str): Position of the variant within the chromosome
        ref (str): The reference allele of the variant
        alt (str): The alternate allele of the variant
    
    Returns:
        list: All variant ids matching the specified variant
    """
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

# CSQ related stuff
def get_csq_field_keys(variantfile):
    """Return an ordered list of csq field keys
    
    Args:
        vcf_filepath (str): Filepath for vcf file to be analysed
    
    Returns:
        List: Consequence field names in the order in which they appear in the vcf
    """

    consequence_tag = "CSQ"

    csq_desc = variantfile.header.info[consequence_tag].description
    csq_field_keys = csq_desc.split("Format: ")[-1].split("|")
    return csq_field_keys

def get_csq_field(field_name, csq, csq_field_keys):
    target_field_index = csq_field_keys.index(field_name)
    csq_fields = csq.split("|")
    csq_target_value = csq_fields[target_field_index]    
    return str(csq_target_value)


def get_csq_dict(variantrecord, csq_field_keys):
    alts = variantrecord.alts
    csqs = variantrecord.info.get("CSQ")
    alleles = alts2csqallele(variantrecord)
    
    csq_dict = {}

    for csq in csqs:
        csq_fields = map(str, csq.split("|"))
        # Which alt does the csq belong to?
        allele_field_name = "Allele"
        allele_field_index = csq_field_keys.index(allele_field_name)
        allele = csq_fields[allele_field_index]
        allele_index = alleles.index(allele)
        alt = alts[allele_index]

        csq_dict.setdefault(alt, []).append(csq)
    
    assert set(csq_dict.keys()) == set(variantrecord.alts), "CSQ missing for allele\n%s " % variantrecord
    #pprint.pprint(csq_dict)
    return csq_dict

    '''
        refseq_field_name = "RefSeq"
        refseq_field_index = csq_field_keys.index(transcript_field_name)
        enst_field_name = "Gene"
        enst_field_index = csq_field_keys.index(enst_field_name)
        
        csq_fields = csq.split("|")
        csq_refseqs = csq_fields[refseq_field_index].split("&")
        csq_ensts = csq_fields[enst_field_index].split("&")

        if csq_transcripts == ['']:
            continue
        print csq_transcripts
        continue
        for gene in gene_data_dict:
            for transcript in gene_data_dict[gene]["transcripts"]:

                # Match when tx versions differ
                if transcript.split(".")[0] in [x.split(".")[0] for x in csq_transcripts]:

                    # We know it's the right transcript
                    # but which alt does the csq belong to
                    allele_field_name = "Allele"
                    allele_field_index = csq_field_keys.index(allele_field_name)
                    allele = csq_fields[allele_field_index]

                    allele_index = alleles.index(allele)
                    alt = alts[allele_index]
                    csq_dict.setdefault(alt, {})
                    csq_dict[alt][transcript] = csq

    #assert len(csq_dict.keys()) > 0, "No relevant csq fields found for variant"
    return csq_dict

    '''

def alts2csqallele(variantrecord):
    """Convert vcf alt field values into equivalent VEP csq allele values
    
    Args:
        vcf_record (_Record): a pyvcf _Record 
    
    Returns:
        TYPE: An ordered list containing one allele for each alt
    """

    csq_alleles = []
    
    for alt in variantrecord.alts:
        
        #alt_seq = alt.sequence

        if len(variantrecord.ref) == 1:
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
                csq_alt = "-"
            else:
            # many:many
                csq_alt = alt[1:]
    
        csq_alleles.append(csq_alt)
    return csq_alleles

# XLS non-variant sheets
def add_summary_sheet(wb, sample_ids, run_folder, styles, gene_data_dict, panels):
    
    ws = wb.add_sheet("Summary")

    # Generate string for panels
    panels_w_id = ["%s (%s)" % (panel_name, pid) 
                for pid, panel_name in panels.items()]
    panels_field = "; ".join(panels_w_id)


    # Calc panel coverage
    total_len = 0
    covered_len = 0

    for gene in gene_data_dict:
        total_len += gene_data_dict[gene]["Total length"]
        covered_len += gene_data_dict[gene]["20x length"]

    coverage_perc = (float(covered_len)/total_len) * 100
    coverage_perc_field = "%0.2f %%" % coverage_perc

    # Infer gender
    sry_field = sry_present(sample, run_folder)
    # TODO - Inferred gender needs to be altered to use het/hom ratio
    if sry_field == "Yes":
        inferred_gender_field = "Male"
    elif sry_field == "No":
        inferred_gender_field = "Female"
    else:
        inferred_gender_field = "Unknown"

    # Report text
    report_text_field = "%s of the target sequence within this panel was sequenced to a depth of 20 fold or more, with analytical sensitivity of 98.3%% - 100%% (95%% Confidence Intervals). The presence of variants reported above, except for variants of unknown significance, has been confirmed by Sanger sequencing. Variants with a population frequency greater than 1 in 500 for dominant conditions, and 1 in 50 for recessive disorders have been deemed insignificant and are not reported. Variants are named using HGVS nomenclature, where nucleotide 1 is the A of the ATG-translation initiation codon." % coverage_perc_field

    # Run QC stats
    read_stats = get_read_stats(sample, run_folder)
    reads_field = "%.2f M." % (read_stats["total"]/1000000.0)
    useable_reads_field = "%.2f M." % ((read_stats["mapped"] - read_stats["duplicate"])/1000000.0)

    # Individual cells

                # Value                 # Style              # Row Col
    cells = [   ("Gemini ID:",          styles["bold"],         0,  0),
                (sample,                styles["plain_text"],   0,  1),
                
                ("GM number:",          styles["bold"],         1,  0),
                #("Value2",             styles["red_text"],     1,  1),
                
                ("Name:",               styles["bold"],         2,  0),
                #("Value3",             styles["red_text"],     2,  1),
            
                ("Inferred Gender:",    styles["bold"],         0,  2),
                (inferred_gender_field, styles["plain_text"],   0,  3),
                
                ("SRY present:",        styles["bold"],         1,  2),
                (sry_field,             styles["plain_text"],   1,  3),

                ("Panel Coverage:",     styles["bold"],         0,  4),
                (coverage_perc_field,   styles["plain_text"],   0,  5),

                ("Panels:",             styles["bold"],         0,  6),
                (panels_field,          styles["plain_text"],   0,  7),

                ("Report Text:",        styles["bold"],         6,  0),
                (report_text_field,     styles["plain_text"],   6,  1),

                ("Bics/Seq QC",         styles["bold"],         35,  0),
                ("Reads:",              styles["bold"],         36,  0),
                (reads_field,           styles["plain_text"],   36,  1),
                ("Usable Reads:",       styles["bold"],         37,  0),
                (useable_reads_field,   styles["plain_text"],   37,  1),

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
    ws = add_table(ws, (27,1), (5,5), column_headers, [(0,3)], styles)
    
    # Column widths
    col_widths = [14,22,22,22,22,22,22,22,12]
    set_col_widths(ws, col_widths)

    return wb

# XLS formatting

def generate_styles():
    # Excel cell formatting - format strings are like css
    style_strings = {"header":"font: bold on;\
                               align: wrap yes;",
                     "bold":"font: bold on",
                     "plain_text":"font: color black",
                     "red_text":"font: color red",
                     "orange_text":"font: color orange",
                     "summary_table_header": "pattern: pattern solid, fore_colour cyan_ega;\
                                              font: color black, bold on;\
                                              borders: left thin, right thin, top thin, bottom thin;",
                     "table_cell":           "borders: left thin, right thin, top thin, bottom thin;",
                     "table_merge_fix":      "borders: left thin;",
                     "hyperlink": "font: color blue",
                    }

    # Use the strings to generate easyxf objects
    styles = {style_name:xlwt.easyxf(style_strings[style_name]) \
                for style_name in style_strings.keys()}

    return styles

def add_table(ws, position, dimensions, col_headers, merged_cols, styles):
    """Add a table to a sheet
    
    Args:
        ws (xlwt.Workbook.Worksheet): an xls worksheet
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

def set_col_widths(ws, col_widths):
    # Set col widths where 1 char ~= 256 units
    for index, col_width in enumerate(col_widths):
        ws.col(index).width = 256 * col_width
    return ws

def set_row_height(ws, row, row_height):
    # Set row height where 1 char ~= 256 units
    ws.row(row).height_mismatch = True  # Unlocks height from cell font height
    ws.row(row).height = 256 * row_height
    return ws

# XLS variant sheet assignment
def generate_variant_type_filter(variant_types):
    def variant_type_filter(csq, csq_field_keys):
        field_name = "Consequence"
        consequences = get_csq_field(field_name, csq, csq_field_keys)
        for variant_type in variant_types:
            if variant_type in consequences:
                return True
        else:
            return False
    return variant_type_filter

def generate_splice_filter(variant_types, flank):
    def splice_filter(csq, csq_field_keys):
        csq_field_name = "Consequence"
        c_dot_field_name = "HGVSc"
        
        consequences = get_csq_field(csq_field_name, csq, csq_field_keys)
        c_dot = get_csq_field(c_dot_field_name, csq, csq_field_keys)

        for variant_type in variant_types:
            if variant_type in consequences:
                # + or - followed by int <= flank
                if int(re.search("[\+-](\d+)", c_dot).groups()[0]) <= flank:
                    return True
        else:
            return False
    return splice_filter

### XLS columns ###
def gene_col(variantrecord, alt, csq, csq_field_keys):
    field_name = "HGNC"
    gene = get_csq_field(field_name, csq, csq_field_keys)
    return gene

def chromosome_col(variantrecord, alt, csq, csq_field_keys):
    return variantrecord.chrom

def start_col(variantrecord, alt, csq, csq_field_keys):
    text = variantrecord.pos
    URL = "https://genome.ucsc.edu/cgi-bin/hgTracks?db=hg19&position=chr%s%%3A%d-%d"\
            % (variantrecord.chrom, variantrecord.pos - 10, variantrecord.pos + 10)
    link = xlwt.Formula('HYPERLINK("%s";"%s")' % (URL, text))

    return text
    return link


def ref_col(variantrecord, alt, csq, csq_field_keys):
    return variantrecord.ref


def alt_col(variantrecord, alt, csq, csq_field_keys):
    return alt


def qual_score_col(variantrecord, alt, csq, csq_field_keys):
    return variantrecord.qual

def dbsnp_col(variantrecord, alt, csq, csq_field_keys):
    
    rsID = variantrecord.id
    return rsID
    
    if rsID:
        ID = rsID[2:]  # remove the rs prefix
        URL = "https://www.ncbi.nlm.nih.gov/projects/SNP/snp_ref.cgi?rs=%s" % ID
        link = xlwt.Formula('HYPERLINK("%s";"%s")' % (URL, rsID))
    
    return link

def refseq_col(variantrecord, alt, csq, csq_field_keys):
    field_name = "RefSeq"
    csq_refseqs = get_csq_field(field_name, csq, csq_field_keys).split("&")
    #assert transcript in csq_transcripts, "Mismatched transcripts!"
    return ", ".join(csq_refseqs)

def enst_col(variantrecord, alt, csq, csq_field_keys):
    field_name = "Gene"
    enst_transcripts = get_csq_field(field_name, csq, csq_field_keys).split("&")
    #assert transcript in csq_transcripts, "Mismatched transcripts!"
    return ", ".join(enst_transcripts)

def c_pos_col(variantrecord, alt, csq, csq_field_keys):
    field_name = "HGVSc"
    tx_c_dot = get_csq_field(field_name, csq, csq_field_keys)
    
    # Sometimes no c dot is present e.g. downstream variants
    if tx_c_dot:
        # Remove tx prefix
        c_dot = tx_c_dot.split(":")[1]
        
        # Find the first letter (excluding c.)
        for index, char in enumerate(c_dot[2:]): # skip the "c." prefix
            if char.isalpha():  
                end_of_pos = index + 2
                break
    
        # Split c_dot into pos and change by splitting at this pos
        c_pos = c_dot[:end_of_pos]
    
    else:
        c_pos = tx_c_dot
    
    return c_pos


def c_change_col(variantrecord, alt, csq, csq_field_keys):
    # Ignored params
    field_name = "HGVSc"
    tx_c_dot = get_csq_field(field_name, csq, csq_field_keys)
    
    # Sometimes no c dot is present e.g. downstream variants
    if tx_c_dot:
        # Remove tx prefix
        c_dot = tx_c_dot.split(":")[1]
        
        # Find the first letter (excluding c.)
        for index, char in enumerate(c_dot[2:]): # skip the "c." prefix
            if char.isalpha():   
                end_of_pos = index + 2
                break
        
        # Split c_dot into pos and change by splitting at this pos
        c_change = c_dot[end_of_pos:]

    else:
        c_change = tx_c_dot
    
    return c_change


def aa_change_col(variantrecord, alt, csq, csq_field_keys):
    # Need to expand aa names from one letter to three letter to be
    # consistent with old report
    field_name = "Amino_acids"
    csq_field = get_csq_field(field_name, csq, csq_field_keys)
    aa_change = csq_field.split("/")
    if len(aa_change) > 1:
        field_name = "Protein_position"
        aa_position = get_csq_field(field_name, csq, csq_field_keys)
        return aa_change[0] + aa_position + aa_change[1]
    else:
        return aa_change[0]

def sift_col(variantrecord, alt, csq, csq_field_keys):
    field_name = "SIFT"
    sift = get_csq_field(field_name, csq, csq_field_keys)
    return sift


def polyphen_col(variantrecord, alt, csq, csq_field_keys):
    field_name = "PolyPhen"
    polyphen = get_csq_field(field_name, csq, csq_field_keys)
    return polyphen

def af_exac_col(variantrecord, alt, csq, csq_field_keys):
    ExAC = "ExAC 0.3"

    vids = get_varfreq_vids(variantrecord.chrom, variantrecord.pos, variantrecord.ref, alt)
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
    return freq
    URL = "http://exac.broadinstitute.org/variant/%s-%s-%s-%s" % (variantrecord.chrom, variantrecord.pos, variantrecord.ref, alt)
    text = freq
    link = xlwt.Formula('HYPERLINK("%s";"%s")' % (URL, text))

    return link

def consequences_col(variantrecord, alt, csq, csq_field_keys):
    field_name = "Consequence"
    consequences = get_csq_field(field_name, csq, csq_field_keys)
    return consequences



@begin.start
def main(vcf_filepath, sample_ids=None, panels=None, maternal=None, paternal=None, index=None, relative=None, exon_flank=50):

    # Either single sample analyses, or a trio
    trio_analysis = all([maternal, paternal, index])
    individual_analysis = not any([maternal, paternal, index, relative])

    assert trio_analysis or individual_analysis, "\nFor trio analysis provide a vcf and maternal, paternal and index.\nFor individual analysis do not use maternal, paternal or index"

    # Check that all input files exist
    filepaths = [f for f in [vcf_filepath] if f]
    for f in filepaths:
        assert os.path.exists(f), "File not found at %s" % f

    run_folder = os.path.basename(os.path.dirname(os.path.dirname(vcf_filepath)))
    output_filepath = get_output_filepath(vcf_filepath)

    # Using pysam to handle the vcf
    vcf_in = pysam.VariantFile(vcf_filepath)
    
    # If no sample ids specified then all in vcf
    # If sample ids specified use those if they are all in vcf
    sample_ids = get_sample_ids(sample_ids, vcf_in)
    
    # Only read the target samples - faster
    vcf_in.subset_samples(sample_ids)

    csq_field_keys = get_csq_field_keys(vcf_in)

    # If no panels specified then use manifest panels
    # If panels specified use those
    sample_panels = sampleids2panels(sample_ids, panels)

    # Sorry future me. This collapses the list of lists of pids into a non-redundant list of pids. 
    # i.e. [[a,b,c],[a,d],[b,e]] -> [a,b,c,d,e]
    unique_pids = list(set([pid for pid_list in [x.keys() for x in sample_panels.values()] for pid in pid_list]))
    
    # Get the genes for each panel
    panel_genes = pids2genes(unique_pids, genepanels_filepath)

    # Collapses the genes in all selected panels into a non-redundant list
    unique_genes = list(set([gene for gene_list in panel_genes.values() for gene in gene_list]))

    # Get the coords (chrom, start, end) for each region associated with
    # the speficied genes in the target bed
    gene_regions = genes2regions(unique_genes, target_bed)

    print "Processing:"
    print "vcf:\t\t%s" % vcf_filepath
    pprint.pprint(sample_panels)
    pprint.pprint(panel_genes)
    pprint.pprint(gene_regions)


    # Build the xls file(s)
    styles = generate_styles()
    wb = xlwt.Workbook()
    #wb = add_summary_sheet(wb, sample_ids, run_folder, styles, gene_data_dict, panel_names)
    if individual_analysis:

        stop_gained_filter = generate_variant_type_filter(["stop_gained_variant","transcript_ablation","start_lost"])
        frameshift_filter =  generate_variant_type_filter(["frameshift_variant","inframe_insertion","inframe_deletion"])
        splice_filter =  generate_splice_filter(["splice_region_variant","splice_acceptor_variant","splice_donor_variant"], 5)
        missense_filter = generate_variant_type_filter(["missense_variant","stop_lost"])
        synonymous_filter = generate_variant_type_filter(["synonymous_variant"])
        # other_filter should always contain "" to match any type and so act
        # as a catch-all for variants which pass through all other filters
        other_filter = generate_variant_type_filter([""])



        sheets = [  #("stop_gained", stop_gained_filter),
                    #("frameshift_variant", frameshift_filter),
                    #("splice", splice_filter),
                    ("missense", missense_filter),
                    #("synonymous", synonymous_filter),
                    ("other", other_filter),
                 ]

        # Specify column headers, functions to generate column values,
        # and col widths (in chars)
        # Cols will appear in the same order in the xls report
        # To add a new column add it to this list

        column_functions =     [("Gene", gene_col, 14),
                                ("RefSeq", refseq_col, 16),
                                ("ENST", enst_col, 16),
                                ("Chr", chromosome_col, 4),
                                ("Start", start_col, 11),
                                ("Ref", ref_col, 5),
                                ("Alt", alt_col, 5),
                                ("Nucleotide pos", c_pos_col, 12),
                                ("Nucleotide change", c_change_col, 12),
                                ("AA change", aa_change_col, 10),
                                ("Score", qual_score_col, 9),
                                #("Depth", depth_col, 7),
                                #("AAF", aaf_col, 6),
                                #("Genotype", genotype_col, 9),
                                ("dbsnp", dbsnp_col, 12),
                                ("Polyphen", polyphen_col, 14),
                                ("SIFT", sift_col, 14),
                                ("Gemini\nAF", af_gemini_col, 9),
                                ("Gemini\nHet/Hom", het_hom_gemini_col, 12),
                                ("ExAC\nAF", af_exac_col, 12),
                                #("ExAC\nHet/Hom", todo_col, 12),
                                #("Gnomad(E)\nAF", af_gnomad_exome_col, 12),
                                #("Gnomad(G)\nAF", af_gnomad_genome_col, 12),
                                ("Consequences", consequences_col, 30),
                                #("OMIM\nIM", OMIM_col, 10), # Use new code for this - see disease forms 
                                ]
    elif trio_analysis:
        sheets = [  ("denovo", denovo_filter),
                    ("recessive", recessive_filter),
                    ("compound_het", comp_het_filter),
                 ]

    else:
        print "OOPS"
        exit()

    sheet_names = [sheet[0] for sheet in sheets]
    column_names = [x[0] for x in column_functions]
    column_widths = [x[2] for x in column_functions]

    # Add variant sheets to the workbook
    variant_sheets = {sheet[0]:wb.add_sheet(sheet[0]) for sheet in sheets}

    # Track the current row in each of the sheets to enable 
    # us to write to the correct line as we write records to 
    # different sheets
    current_row_index = {sheet[0]:0 for sheet in sheets}
    
    # Add variant sheet headers
    header_style = styles["header"]
    for sheet_name, ws in variant_sheets.items():
        for col_index, col in enumerate(column_functions):
            col_name, column_function, _ = col
            row_index = current_row_index[sheet_name]
            ws.write(row_index, col_index, col_name, style=header_style)
        
        freeze_header(ws)  # Keep header on screen when scrolling
        set_col_widths(ws, column_widths)
        header_row_height = 3
        set_row_height(ws, 0, header_row_height)
        current_row_index[sheet_name] += 1

    wb.save("test.xls")

    # Identify the variants affecting the gene
    for gene, gene_data in sorted(gene_regions.items()):
        #print gene
        gene_variants = []

        for transcript, transcript_data in sorted(gene_data.items()):
            #print transcript
            for exon, exon_data in sorted(transcript_data.items()):
                #print exon
                for variant in vcf_in.fetch(exon_data["chrom"], exon_data["start"]-exon_flank, exon_data["end"]+exon_flank):
                    # These are our variants of interest
                    #print variant
                    gene_variants.append(variant)

        print
        print gene
        #print gene_variants

        # Might need to set() gene_variants here

        
        for v_index, variant in enumerate(gene_variants):
            print variant
            print "V-A-C"
            var_prefix = str(v_index)
            csq_dict = get_csq_dict(variant, csq_field_keys)    
            
            for a_index, alt in enumerate(variant.alts):
                csqs = csq_dict[alt]
                alt_prefix = "-"+str(a_index)
                
                for c_index, csq in enumerate(csqs):
                    csq_prefix = "-"+str(c_index)
                    row_data = []

                    for c in column_functions:
                        column_header = c[0]
                        column_function = c[1]
                        column_name, column_value = \
                            column_header,\
                            column_function(variantrecord=variant,
                                            alt=alt, 
                                            #transcript=transcript, 
                                            csq=csq,
                                            csq_field_keys=csq_field_keys,
                                            )
                        row_data.append(column_value)
                    print var_prefix+alt_prefix+csq_prefix, row_data
                #print list(variant.info)

