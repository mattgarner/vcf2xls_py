"""Summary
"""
import os
import sys
import MySQLdb
import begin
import vcf
import xlwt
import urllib2
import subprocess
import pprint


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


def sry_depth(sample):
    query = """
            SELECT mean
            FROM coverage
            WHERE sid = (SELECT sid 
                         FROM sample 
                         WHERE name = '%s')
            AND rid = (SELECT rid
                       FROM region2gene
                       WHERE gid = (SELECT gid
                                    FROM gene
                                    WHERE name = "SRY"))
            """ % sample

    depth = geminiDB_query(query)[0][0]
    return depth


def get_read_stats(sample):
    query = """
            SELECT total_reads, mapped_reads, duplicate_reads
            FROM stats
            WHERE sid = (SELECT sid
                         FROM sample
                         WHERE name = '%s')
            """ % sample

    read_stats = geminiDB_query(query)[0]

    total, mapped, duplicate = read_stats
    read_stats = {}
    read_stats["total"] = total
    read_stats["mapped"] = mapped
    read_stats["duplicate"] = duplicate

    return read_stats


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
    pid2panel = {}
    
    for pid in pids:
        panel_name = pid2panelname(pid, genepanels_filepath)
        pid2panel[pid] = panel_name

    return pid2panel


def pid2panelname(pid, genepanels_filepath):
    if pid.startswith("_"):
        panel_name = pid
    elif pid == "8926":
        panel_name = "TSO"
    
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


def panelnames2pids(panel_names, genepanels_filepath):
    panel2pid = {}
    for panel_name in panel_names.split(","):
        pid = panelname2pid(panel_name, genepanels_filepath)
        panel2pid[panel_name] = pid
    return panel2pid


def panelname2pid(panel_name, genepanels_filepath):
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


def pids2genes(pids, genepanels_filepath, genes2transcripts_filepath):
    """Return a dict of panel_ids and their genes
    
    Args:
        panels (list): A list of panel_IDs
        genepanels_filepath (str): Filepath for tab delimited text file containing Panel_Name Panel_ID Gene_Name
    
    Returns:
        TYPE: A alphabetically sorted non-redundant list of names of all genes in all specified panels
    """
    
    panel_genes = {}

    # Some genes in the genepanels file have old HGNC identifiers. 
    # Consequently these gene names do not match those in the 
    # genes2transcripts file.
    # This dict is used to map the old alias (the key) to the current
    # id in the genes2transcripts file(the value).

    gene_alias_dict = get_gene_alias_dict()

    for pid in pids:
        print "pid", pid
        genes = pid2genes(pid, genepanels_filepath, genes2transcripts_filepath)
        # Replace old aliases
        genes = [gene_alias_dict.get(gene, gene) for gene in genes]
        print genes
        panel_genes[pid] = genes
    
    print panel_genes
    return panel_genes


def pid2genes(pid, genepanels_filepath, genes2transcripts_filepath):
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
            for line in genelist_fh.readlines()[:40]:  # DEBUG
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


def genes2transcripts(genes, genes2transcripts_filepath):
    """Return a dict of gene_name:[transcripts] for all genes
    
    Args:
        genes (list): A list of HGNC gene names
        genes2transcripts_filepath (str): Filepath for tab delimited text file containing Gene_Name Transcript
    
    Returns:
        dict: gene_name:[transcripts]
    """

    gene_transcripts = {}
    
    for gene in genes:
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
    
    return list(set(transcripts_for_gene))


def transcripts2exons(transcripts):
    transcript_exons = {}
    
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


def record_in_regions(vcf_record, gene_data_dict, flank=0):
    vcf_record_range = record_range(vcf_record)
    
    region_hits = []

    for gene in gene_data_dict:
        for transcript in gene_data_dict[gene]["transcripts"]:

            for exon, exon_data in sorted(gene_data_dict[gene]["transcripts"][transcript]["exons"].items()):

                # Wrong chrom
                if vcf_record.CHROM != exon_data["chr"]:
                    break
                
                # Right chrom - do coords overlap?
                else:
                    region_range = [exon_data["start"]-flank, exon_data["end"]+flank]
                    
                    overlap_length = overlap(region_range, vcf_record_range)
                    if overlap_length:
                        #print "v", vcf_record_range
                        #print "p", region_range
                        #print overlap_length
                        region_hits.append(transcript)
                    
    return region_hits
                   

def record_range(vcf_record):
    alleles = vcf_record.ALT + [vcf_record.REF]
    max_allele_length = max([len(allele) for allele in alleles])
    
    first = vcf_record.POS
    last = vcf_record.POS + (max_allele_length -1)

    return [first, last]


def overlap(a, b):
    a = map(int, [x for x in a])
    b = map(int, [x for x in b])

    assert all([x >= 0 for x in a]), "Range cannot contain negative number"
    assert all([x >= 0 for x in b]), "Range cannot contain negative number"
    return max(0, min(a[1]+1, b[1]+1) - max(a[0], b[0]))


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


def get_gene_alias_dict():
    gene_alias_dict = {   'DTDST'      :'SLC26A2',
                          'EKLF'       :'KLF1',
                          # chromosome breaking sydromes
                          'C7orf11'    :'MPLKIP',
                          'C7ORF11'    :'MPLKIP',
                          'CSA/ERCC8'  :'ERCC8',
                          'TTDA/GTF2H5':'GTF2H5',
                          'XPG/ERCC5'  :'ERCC5',
                          'DDB2/XPE'   :'DDB2',
                          'CSB/ERCC6'  :'ERCC6',
                          'XPB/ERCC3'  :'ERCC3',
                          'XPD/ERCC2'  :'ERCC2',
                          'XPF/ERCC4'  :'ERCC4',
                          'XPV/POLH'   :'POLH',
                          'FANCG/XRCC9':'FANCG',
                          'FANCL/PHF9' :'FANCL',
                          'OX40'       :'TNFRSF4',
                          'MRE11'      :'MRE11A',
                          'IGHM'       :'IGHMBP2',
                          'FATP4 (SLC27A4)':'SLC27A4',
                          'FATP4(SLC27A4)' :'SLC27A4',
                          'FATP4'      :'SLC27A4',
                          'KCNE5'      :'KCNE1L',
                          'P3H1'       :'LEPRE1',
                          'GNAS1'      :'GNAS',
                          'FBLN4'      :'EFEMP2',
                          }
    return gene_alias_dict



#### Report columns ####

# functions should be structured like:
# def field_col(vcf_record, alt, transcript, csq, csq_field_keys)
    # Ignored params
    # del ...
    # do stuff
    # return value_for_col


# vcf field functions
# These only require a vcf_record as input
# Generally DNA level info
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


# csq field functions
# There require a specific vcf_record csq as input
# Generally RNA/protein level info
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


def c_change_col(vcf_record, alt, transcript, csq, csq_field_keys):
    # Ignored params
    del vcf_record, alt, transcript
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

    genemap = get_genemap(gene)

    if not genemap:
        return ""

    mim = genemap["geneMap"]["mimNumber"]

    phenotype_inheritance = get_phenotype_inheritance(genemap)

    URL = "https://www.omim.org/entry/%s" % mim
    text = phenotype_inheritance
    link = xlwt.Formula('HYPERLINK("%s";"%s")' % (URL, text))

    return link


# Other field functions
# These require the specific alt and/or transcript as input
def transcript_col(vcf_record, alt, transcript, csq, csq_field_keys):
    # Ignored params
    del vcf_record, alt
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


def todo_col(vcf_record, alt, transcript, csq, csq_field_keys):
    return "TODO"

# Non-variant sheets
def add_summary_sheet(wb, sample, styles, gene_data_dict, panels):
    
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
    sry = sry_present(sample)
    if sry:
        sry_field = "Yes"
        inferred_gender_field = "Male"
    else:
        sry_field = "No"
        inferred_gender_field = "Female"

    # Report text
    report_text_field = "%s of the target sequence within this panel was sequenced to a depth of 20 fold or more, with analytical sensitivity of 98.3%% - 100%% (95%% Confidence Intervals). The presence of variants reported above, except for variants of unknown significance, has been confirmed by Sanger sequencing. Variants with a population frequency greater than 1 in 500 for dominant conditions, and 1 in 50 for recessive disorders have been deemed insignificant and are not reported. Variants are named using HGVS nomenclature, where nucleotide 1 is the A of the ATG-translation initiation codon." % coverage_perc_field

    # Run QC stats
    read_stats = get_read_stats(sample)
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


def add_geneqc_sheet(wb, sample, styles, gene_data_dict):
    
    ws = wb.add_sheet("Gene_QC")

    columns = [ "Gene","Transcript","Region length","Min depth",
                "Missing","1-5x","6-9x","10-19x","20+x","20+x %" ]
    
    row_index = 0

    # Write headers
    for col_index, column_header in enumerate(columns):
        ws.write(row_index, col_index, column_header, style=styles["bold"])
    row_index += 1

    total_len = sum([gene_data["Total length"] for gene_data in gene_data_dict.values()])
    twentyx_len = sum([gene_data["20x length"] for gene_data in gene_data_dict.values()])
    twentyx_perc = (float(twentyx_len)/total_len) * 100

    for gene in sorted(gene_data_dict):
        for transcript in gene_data_dict[gene]["transcripts"]:

            # Set style rules here
            # Default
            row_style = styles["plain_text"]
            # Exceptions
            if float(gene_data_dict[gene]["transcripts"][transcript]["20+x %"]) < 100.0:
                row_style = styles["red_text"]


            # Write data rows
            for col_index, column_header in enumerate(columns):
                value = gene_data_dict[gene]["transcripts"][transcript][column_header]
                ws.write(row_index, col_index, value, row_style)
            row_index += 1
        row_index += 1

    # The 'Total' row
    ws.write(row_index, 0, "Total:")
    ws.write(row_index, columns.index("Region length"), total_len)
    ws.write(row_index, columns.index("20+x"), twentyx_len)
    ws.write(row_index, columns.index("20+x %"), "%.2f" % twentyx_perc)

    ws = freeze_header(ws)

    # Column widths
    col_widths = [14,16,14,10,8,8,8,8,8,8,8,8]
    set_col_widths(ws, col_widths)

    # Header row height
    set_row_height(ws, 0, 3)

    return wb


def add_exonqc_sheet(wb, sample, styles, gene_data_dict):
    
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
    for gene in sorted(gene_data_dict):

        for transcript in gene_data_dict[gene]["transcripts"]:

            # Write data rows
            for exon, exon_data in sorted(gene_data_dict[gene]["transcripts"][transcript]["exons"].items()):

                # Set style rules here
                # Default
                row_style = styles["plain_text"]
                # Exceptions
                if int(exon_data["Min depth"]) < 20:
                    row_style = styles["red_text"]
                if exon_data["CNV type"]:
                    row_style = styles["orange_text"]
                
                for col_index, column_header in enumerate(columns):
                    value = exon_data.get(column_header, "")
                    ws.write(row_index, col_index, value, row_style)
                row_index += 1
            row_index += 1

    ws = freeze_header(ws)

    # Column widths
    col_widths = [14,16,30,12,12,12,12,8,8,8,8,12,24,10,16]
    set_col_widths(ws, col_widths)

    # Header row height
    set_row_height(ws, 0, 3)

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


# Sheet formatting
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


def freeze_header(ws):
    ws.set_panes_frozen(True)
    ws.set_horz_split_pos(1)
    return ws


def set_col_widths(ws, col_widths):
    # Set col widths - 1 char ~= 256 units
    for index, col_width in enumerate(col_widths):
        ws.col(index).width = 256 * col_width
    return ws


def set_row_height(ws, row, row_height):
    # Set row height - 1 char ~= 256 units
    ws.row(row).height_mismatch = True  # Unlocks height from cell font height
    ws.row(row).height = 256 * row_height
    return ws


def write_row_to_worksheet(ws, row_number, row_data, style):
    for col_index, value in enumerate(row_data):
        # print "Row:\t%d\tCol:\t%d\tValue:\t%s" % (row_number, col_index, value)
        ws.write(row_number, col_index, value, style)


# CSQ stuff
def get_csq_field(field_name, csq, csq_field_keys):

    target_field_index = csq_field_keys.index(field_name)
    csq_fields = csq.split("|")
    csq_target_value = csq_fields[target_field_index]    
    return csq_target_value


def get_csq_dict(vcf_record, gene_data_dict, csq_field_keys):
    
    alts = vcf_record.ALT
    csqs = vcf_record.INFO.get("CSQ")
    alleles = alts2csqallele(vcf_record)

    csq_dict = {}

    for csq in csqs:
        transcript_field_name = "RefSeq"
        transcript_field_index = csq_field_keys.index(transcript_field_name)
        
        csq_fields = csq.split("|")
        csq_transcripts = csq_fields[transcript_field_index]
        
        for gene in gene_data_dict:
            for transcript in gene_data_dict[gene]["transcripts"]:
                if transcript in csq_transcripts:
                    print
                    print vcf_record
                    print csq
                    print alleles
                    print gene, transcript
                    # it's the right transcript
                    # which alt does it belong to
                    allele_field_name = "Allele"
                    allele_field_index = csq_field_keys.index(allele_field_name)
                    allele = csq_fields[allele_field_index]
                    print allele
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
                csq_alt = alt_seq[1:]
    
        csq_alleles.append(csq_alt)
    return csq_alleles


# Decon for dosage
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
def get_gene_data_dict(sample, panel_ids, genepanels_filepath, genes2transcripts_filepath, decon_filepath):
    
    print 1
    panel_names = pids2panelnames(panel_ids, genepanels_filepath)
    print 2
    panel_genes = pids2genes(panel_ids, genepanels_filepath, genes2transcripts_filepath)
    print 3
    gene_names = sorted([gene for genes_list in panel_genes.values()
                            for gene in genes_list])
    print 4
    gene_transcripts = genes2transcripts(gene_names, genes2transcripts_filepath)
    print 5 
    transcript_ids = [transcript for transcripts_list in gene_transcripts.values()
                            for transcript in transcripts_list]
    print 6
    transcript_regions = transcripts2exons(transcript_ids)
    print 7

    #decon_filepath = get_decon_results_filepath(sample)
    decon_data_dict = decon_dict(decon_filepath)
    print 8

    gene_data_dict = {}

    # Define genes, transcripts, exons in panel
    
    for panel_id in panel_ids:
        print "PROCESSING PANEL:", panel_id
        
        for gene in panel_genes[panel_id]:
            print gene

            gene_data_dict.setdefault(gene, {})
            
            for transcript in gene_transcripts[gene]:

                gene_data_dict[gene].setdefault("transcripts", {})
                gene_data_dict[gene]["transcripts"][transcript] = {}
                                
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


                #first_record = True
                for record in result:
                    gene_data_dict[gene]["transcripts"][transcript].setdefault("exons", {})
                    
                    # We get one record per exon
                    sample, HGNC, exon, refseq, rid, chr_, start, end, min_, mean,\
                    max_, missing, oneto5, sixto9, tento19 = record

                    # Pull values straight from query
                    gene_data_dict[gene]["transcripts"][transcript]["exons"][exon] = \
                        {  "Sample":sample,
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
                    gene_data_dict[gene]["transcripts"][transcript]["exons"][exon]["Name"] = calc_exon_name(gene_data_dict[gene]["transcripts"][transcript]["exons"][exon])
                    gene_data_dict[gene]["transcripts"][transcript]["exons"][exon]["Position"] = calc_exon_position(gene_data_dict[gene]["transcripts"][transcript]["exons"][exon])
                    gene_data_dict[gene]["transcripts"][transcript]["exons"][exon]["Exp mean depth"] = calc_exon_exp_mean_depth(gene_data_dict[gene]["transcripts"][transcript]["exons"][exon])

                    # Get DECoN dosage calls for exons - adds multiple fields
                    gene_data_dict[gene]["transcripts"][transcript]["exons"][exon] = update_exon_with_decon(gene_data_dict[gene]["transcripts"][transcript]["exons"][exon], decon_data_dict)


    # Once the full exon dict is built we then generate values for each
    # gene(# missing bases etc) by using the dict

    # The order of these calls matters - some are dependent on the
    # output of preceeding calls being present in the data dict
    print "\n\n\nAdding gene stats"
    for gene in sorted(gene_data_dict):
        print gene
        transcripts = gene_data_dict[gene].get("transcripts")
        assert transcripts, "No Tx for %s. Aborting!" % gene
        for transcript in transcripts:
            print transcript

            transcript_data = gene_data_dict[gene]["transcripts"][transcript]
            
            transcript_data["Gene"] = gene
            transcript_data["Transcript"] = transcript
            transcript_data["Region length"] = calc_region_length(transcript_data)
            transcript_data["Min depth"] = calc_gene_min_depth(transcript_data)
            transcript_data["Missing"] = calc_bases_in_depth_bin(transcript_data, "0x")
            transcript_data["1-5x"] = calc_bases_in_depth_bin(transcript_data, "1-5x")
            transcript_data["6-9x"] = calc_bases_in_depth_bin(transcript_data, "6-9x")
            transcript_data["10-19x"] = calc_bases_in_depth_bin(transcript_data, "10-19x")
            transcript_data["20+x"] = calc_gene_20x_len(transcript_data)
            transcript_data["20+x %"] = calc_gene_20x_perc(transcript_data)    
    
    # If we used one tx per gene then this section would
    # not be necessary - we could just use the tx stats.
    # Since we use multiple txs in some cases we need 
    # stats for the region comprising all collapsed
    # transcripts we use for each gene
    
        # For each gene
        # x = what is the total length (merge all exons intervals)
        # y = what is the total non-20x (merge all sub20x intervals)
        # y/x * 100 = 20x % for gene
        gene_data = gene_data_dict[gene]
        exon_ranges = []
        sub_20x_ranges = []

        for transcript, transcript_data in gene_data["transcripts"].items():
            for exon, exon_data in transcript_data["exons"].items():
                exon_ranges.append((exon_data["start"], exon_data["end"]))
        
                # Each exon data field can contain more than one range
                # e.g. "1:123-125, 1:130-132"
                # Split these into individual ranges on ","
                # Then merge all the split ranges into one list
                # Makes a list like [1:123-456, 1:789-1011]
                # First merge the strings into one.Then split on "," to 
                # separate ranges into a list with one range per item
        
                lt_20x = ",".join([exon_data["0x"    ],
                                   exon_data["1-5x"  ],
                                   exon_data["6-9x"  ],
                                   exon_data["10-19x"]
                                   ])

                lt_20x = [x for x in lt_20x.split(",") if x != ""]
                
                for lt_20x_range in lt_20x:
                    start, end = map(int, lt_20x_range.split(":")[1].split("-"))
                    sub_20x_ranges.append((start, end))

        merged_sub_20x_ranges = merge_intervals(sub_20x_ranges)
        sub_20x_len = sum([(r[1]-r[0])+1 for r in merged_sub_20x_ranges])

        merged_exon_ranges = merge_intervals(exon_ranges)
        total_len = sum([(r[1]-r[0])+1 for r in merged_exon_ranges])

        gene_data = gene_data_dict[gene]
        
        gene_data["Total length"] = total_len
        gene_data["20x length"] = total_len - sub_20x_len
        gene_data["20+x %"] = (float(gene_data["20x length"])/total_len) * 100

    '''
    gene_data = gene_data_dict[gene]
    gene_data["Gene"] = gene
    gene_data["Transcript"] = gene_data_dict[gene]["transcripts"].keys()
    gene_data["Collapsed regions"] = []
    gene_data["Region length"] = calc_region_length(transcript_data)
    #gene_data["Min depth"] = calc_gene_min_depth(transcript_data)
    #gene_data["Missing"] = calc_bases_in_depth_bin(transcript_data, "0x")
    #gene_data["1-5x"] = calc_bases_in_depth_bin(transcript_data, "1-5x")
    #gene_data["6-9x"] = calc_bases_in_depth_bin(transcript_data, "6-9x")
    #gene_data["10-19x"] = calc_bases_in_depth_bin(transcript_data, "10-19x")
    gene_data["20+x"] = calc_gene_20x_len(transcript_data)
    gene_data["20+x %"] = calc_gene_20x_perc(transcript_data)                                    
    '''
    
    return gene_data_dict
           

def calc_region_length(transcript_data_dict):
    region_length = sum([exon["end"] - exon["start"] + 1 for exon in transcript_data_dict["exons"].values()])
    return region_length


def calc_gene_min_depth(transcript_data_dict):
    gene_min_depth = float("inf")

    for exon in transcript_data_dict["exons"].values():
        exon_min_depth = exon["Min depth"]
        if (gene_min_depth == None) or (gene_min_depth > exon_min_depth):
            gene_min_depth = exon_min_depth

    if gene_min_depth == float("inf"):
        gene_min_depth = None 

    return gene_min_depth


def calc_bases_in_depth_bin(transcript_data_dict, depth_bin):
    
    bases_in_depth_bin = 0

    for exon in transcript_data_dict["exons"].values():
        regions = exon[depth_bin].split(",")
        for region in regions:
            if region:
            # e.g. 7:117232563-117232564
                start, end = map(int, region.split(":")[1].split("-"))
                bases_in_depth_bin += end-start+1

    return bases_in_depth_bin


def calc_gene_20x_len(transcript_data_dict):
    # No 20x+ bin in database, so need to calc this using other bins

    total_len = transcript_data_dict["Region length"]
    under_20x_len = 0

    for exon in transcript_data_dict["exons"].values():
        for depth in ["0x", "1-5x", "6-9x", "10-19x"]:
            low_depth_regions = exon[depth].split(",")
            for region in low_depth_regions:
                if region:
                # e.g. 7:117232563-117232564
                    start, end = map(int, region.split(":")[1].split("-"))
                    under_20x_len += end-start+1
    len_20x = total_len - under_20x_len

    return len_20x


def calc_gene_20x_perc(transcript_data_dict):
   len_total = transcript_data_dict["Region length"]
   len_20x = transcript_data_dict["20+x"]
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


def merge_intervals(ranges):
    # create a list of starts and ends
    stends = [ (r[0], 1) for r in ranges ]
    stends += [ (r[1], -1) for r in ranges ]
    stends.sort(key=lambda x: x[0])

    # now we should have a list of starts and ends ordered by position,
    # e.g. if the ranges are 5..10, 8..15, and 12..13, we have
    # (5,1), (8,1), (10,-1), (12,1), (13,-1), (15,-1)

    # next, we form a cumulative sum of this
    s = 0
    cs = []
    for se in stends:
        s += se[1]
        cs.append((se[0], s))

    # this is, with the numbers above, (5,1), (8,2), (10,1), (12,2), (13,1), (15,0)
    # so, 5..8 belongs to one range, 8..10 belongs to two overlapping range,
    # 10..12 belongs to one range, etc

    # now we'll find all contiguous ranges
    # when we traverse through the list of depths (number of overlapping ranges), a new
    # range starts when the earlier number of overlapping ranges has been 0
    # a range ends when the new number of overlapping ranges is zero 
    prevdepth = 0
    start = 0
    merged = []
    for pos, depth in cs:
        if prevdepth == 0:
            start = pos
        elif depth == 0:
            merged.append((start, pos))
        prevdepth = depth

    return merged


# Variant sheet filters
def generate_variant_type_filter(variant_types):
    def variant_type_filter(csq, csq_field_keys):
        field_name = "Consequence"
        consequences = get_csq_field(field_name, csq, csq_field_keys)
        #print consequences
        for variant_type in variant_types:
            #print "Searching for:\n%s\nin\n%s\n" % (variant_type, consequences)
            if variant_type in consequences:
            #    print True
                return True
        else:
            #print False
            return False
    return variant_type_filter

stop_gained_filter = generate_variant_type_filter(["stop_gained_variant","transcript_ablation","start_lost"])
frameshift_filter =  generate_variant_type_filter(["frameshift_variant","inframe_insertion","inframe_deletion"])
splice_filter =  generate_variant_type_filter(["splice_region_variant","splice_acceptor_variant","splice_donor_variant"])
missense_filter = generate_variant_type_filter(["missense_variant","stop_lost"])
synonymous_filter = generate_variant_type_filter(["synonymous_variant"])
# other_filter should always contain "" to act as a catch-all for variants
# which pass all other filters
other_filter = generate_variant_type_filter([""])


#### OMIM Stuff ####
"""
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
"""

def construct_api_request(search_term):
    '''
    Assemble an OMIM API URL for a specific search term
    '''
    omim_api_url = "http://api.omim.org/api"
    response_format = "python"
    handler = "geneMap"
    action = "search"
    API_key = "apiKey=YwfXGOC4Sim_82R9qt8Xzw"
    limit_results = 5
    search_params = "search=gene_symbol%3B" + str(search_term)
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

"""
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

"""

def get_genemap(query_HUGO):
    '''
    Extract geneMap for a specific HGNC gene name from returned OMIM API data
    '''
    query_HUGO = query_HUGO.upper()
    print query_HUGO

    genemaps = []

    api_request = construct_api_request(query_HUGO)
    api_response = submit_api_request(api_request)

    # Grab the genemap specific for the gene - search may return others
    search_results = api_response["omim"]["searchResponse"]["geneMapList"]
    
    for genemap in search_results:
        
        # OMIM symbols are mixed case, 
        gene_symbols = [x.upper() for x in genemap["geneMap"]["geneSymbols"].split(",")]

        if query_HUGO in gene_symbols:
            genemaps.append(genemap)

    #assert len(genemaps) <= 1, "UNEXPECTED NUMBER OF GENEMAPS %d" % len(genemaps)

    if len(genemaps) != 1:
        return None
    return genemaps[0]


def get_phenotype_inheritance(genemap):
    if genemap == None:
        return ""

    phenotype_inheritance = {}
 
    # Phenotype level data - many per gene/entry
    for phenotype in genemap["geneMap"].get("phenotypeMapList", []):
        phenotype_name = phenotype.get("phenotypeMap", {}).get("phenotype")
        phenotype_inheritance_mode = phenotype.get("phenotypeMap", {}).get("phenotypeInheritance")
        
        # Don't return 'nondiseases' ' [ ' 
        # or susceptibility for multifactorial diseases ' { 
        if phenotype_name[0] not in ["[", "{"]:
            if phenotype_inheritance_mode != None:
                phenotype_inheritance[phenotype_name] = phenotype_inheritance_mode

    # Collapse results
    all_pheno_inheritance = []
    for phenotype in phenotype_inheritance:
        inheritance_modes = list(set(x for x in phenotype_inheritance.values()))
        all_pheno_inheritance += inheritance_modes
    collapsed_pheno_inheritance = list(set(all_pheno_inheritance))

    # Abbreviate terms
    inheritance_mode_mapping = [("Autosomal recessive","AR"),
                                ("Autosomal dominant","AD"),
                                ("X-linked dominant","XD"),
                                ("X-linked recessive","XR"),
                                ]

    for index, inheritance_mode in enumerate(collapsed_pheno_inheritance):
        for long_name, short_name in inheritance_mode_mapping:
            collapsed_pheno_inheritance[index] = \
                collapsed_pheno_inheritance[index].replace(long_name, short_name)
    
    return ",".join(collapsed_pheno_inheritance)


# Other
def get_vcf_sample_ids(vcf_filepath):
    """Return a list of all sample ids in a vcf file
    
    Args:
        vcf_filepath (str): Filepath to a vcf file to be evaluated
    
    Returns:
        list: Sample ids
    """
    with open(vcf_filepath) as vcf_fh:
        vcf_reader = vcf.Reader(vcf_fh)
        sample_ids = vcf_reader.samples
    return sample_ids


def sry_present(sample, depth_threshold=50):
    sry_read_depth = sry_depth(sample)
    if sry_read_depth >= depth_threshold:
        return True
    else:
        return False


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


def annotate_vcf(vcf_filepath):
    if vcf_filepath.endswith(".annotated.vcf"):
        return vcf_filepath
    else:
        print "TODO - NEED TO ANNOTATE UNANNOTATED VCFS"
        print "USE NEW PIPELINE?"
        exit()


def gene_coverage(gene_data_dict):
    pass


@begin.start
def main(
    vcf_filepath,
    decon_filepath="/data/projects/matt/cnv/DECoN_results/170405_170407_NS500192/170405_170407_NS500192_all.txt",
    samplepanels_filepath="/data/gemini/BioinformaticManifest.txt",
    genes2transcripts_filepath="/data/gemini/genes2transcripts",
    genepanels_filepath="/data/gemini/genepanels",
    panels=None,
    exon_flank=30
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

    # Check that the vcf file contains data for the sample specified 
    # in the filename
    vcf_sample_ids = get_vcf_sample_ids(vcf_filepath)
    assert len(vcf_sample_ids) == 1,\
        "Multiple samples in vcf file"

    sample_id = os.path.basename(vcf_filepath).split(".")[0]

    assert sample_id in vcf_sample_ids,\
        "SampleID in filename not found in file!"

    output_filepath = get_output_filepath(vcf_filepath)

    # Annotate any unannotated vcfs
    vcf_filepath = annotate_vcf(vcf_filepath)


    #### Generate sample/panel specific data ####
    print "### Generating panel data ###"
    if panels:
        if panels == "full":
            panel_IDs = ["8926"]
        else:
            panel_IDs = panelnames2pids(panels, genepanels_filepath)
    else:
        panel_IDs = gnum2pids(sample_id, samplepanels_filepath)

    gene_data_dict = get_gene_data_dict(sample_id, panel_IDs, genepanels_filepath, genes2transcripts_filepath, decon_filepath)

    panel_names = pids2panelnames(panel_IDs, genepanels_filepath)


    #### Generate vcf specific metadata ####
    csq_field_keys = get_csq_field_keys(vcf_filepath)


    # Set up the xls workbook
    with open(vcf_filepath) as vcf_fh:
        # strict_whitespace=True splits on tabs only, enabling spaces
        # to be present in fields which shouldn't really have them 
        # according to vcf spec.
        # Our vcfs do have such spaces(in the VEP consequence strings)
        # so this is required to correctly parse the vcfs
        vcf_reader = vcf.Reader(vcf_fh, strict_whitespace=True)

        styles = generate_styles()

        # Make a workbook
        wb = xlwt.Workbook()

        # Add the non-variant sheets to the workbook
        print "### Adding QC sheets ###"
        wb = add_summary_sheet(wb, sample_id, styles, gene_data_dict, panel_names)
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
        # To add a new column add it to this list, along with a *_col
        # function to generate the value

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
                                #("OMIM\nIM", OMIM_col, 10),
                                ]

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
        
        # Tracks the gene of the most recently added record
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
                            current_row_index[ws.name] += 1
                            current_gene_per_sheet[ws.name] = row_data[column_names.index("Gene")]
                            break

                    
    print "### Report complete! ###"
    wb.save(output_filepath)
    print "Output: %s" % output_filepath

