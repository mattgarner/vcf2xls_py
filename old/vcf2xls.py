import os
import sys
import xlwt
import begin
import vcf


def get_sampleID(sample):
    pass


def get_run_folder(sample):
    pass


def get_gm_number(sample):
    pass


def get_patient(sample):
    pass


def open_fh(filepath, read_write="r", ext=None):
    """More robust filehandle opening
    
    Args:
        filepath (string): Filepath of the file to be opened
        read_write (str, optional): Specify whether fh is for reading or writing the file.
        ext (None, optional): Specify an particular extension for the filepath
    
    Returns:
        file: A filehandle for the specified file, in the specified mode
    """

    # Filepath exists
    assert os.path.exists(filepath),\
        "File not found: %s" % filepath
    
    # Filepath is a file
    assert os.path.isfile(filepath),\
        "Filepath must point to a file, not a directory: %s" % filepath
    
    # File has correct ext (if specified)
    if ext:
        # In case the . was not included in the specified ext, add it
        if not ext.startswith("."):
            ext = "." + ext
        assert os.path.splitext(filepath)[1] == ext,\
            "File %s does not have extension: %s" % (filepath, ext)
    
    # User has specified a valid mode for opening the file
    assert read_write in ["r","w"], "read_write must be one of 'r'(ead) or 'w'(rite)"
    
    # Everything looks good, open the filehandle
    filehandle = open(filepath, read_write)
    return filehandle


def is_variant_type(line):
    pass


def generate_variant_sheets(vcf_filepath):
    report = xlwt.Workbook()
    variant_types = ["stop_gained","frameshift_variant","consensus_splice","missense_variant","synonymous", "other"]
    # might need to make each a dict of key = name, value = list of seq ontology terms for type
    for variant_type in variant_types[0:1]:
        ws = wb.add_sheet(variant_type)
        generate_variant_sheet(variant_type, vcf_filepath, ws)


def generate_variant_sheet(variant_type, vcf_filepath, worksheet):
    vcf_reader = vcf.Reader(open(vcf_filepath, 'r'))
    csq_fields = get_CSQ_format(vcf_reader)
    sheet_columns = 
    # Maybe each col maps to a func which pulls the value for that col
    # Adding new col just requires adding a new col, a new field in vcf, and a function mapping the two
    # Makes future dev much more straighforward
    # Probably a bit slower than a crazy script

    for row_index, vcf_record in enumerate(vcf_reader):
        csq_list = vcf_record.INFO.get("CSQ")
        for csq_string in csq_list:
            csq_dict = csq2dict(csq_string, csq_fields)
            print vcf_record
            print csq_dict
            #raw_input()


def get_CSQ_format(vcf_reader):
    format_tag = "Format: "
    csq_desc = vcf_reader.infos.get("CSQ").desc
    csq_format_pos = csq_desc.find(format_tag)
    csq_fields = csq_desc[csq_format_pos+len(format_tag):].split("|")
    return csq_fields


def csq2dict(csq_string, csq_fields):
    csq_dict = {}
    csq_values = csq_string.split("|")
    
    # List of tuples, each containing a key and it's value
    csq_key_value_pairs = zip(csq_fields, csq_values)
    
    # Make a dict of the pairs
    for pair in csq_key_value_pairs:
        key, value = pair
        csq_dict[key] = value

    return csq_dict


'''
def generate_variant_sheet(variant_type, vcf_filepath):
    vcf_fh = open_fh(vcf_filepath, "w", ".vcf")
    vcf_fileheader = []
    
    for line in vcf_fh:
        line = line.strip()
        if line.startswith("##"):
            vcf_fileheader.append(line)
        elif line.startswith('#'):
            vcf_col_headers = line.split()
        else:
            vcf_line = parse_vcf_line(line, vcf_col_headers)
            if is_variant_type(vcf_line, variant_type):
            process_vcf_line(line, )
'''

@begin.start
def main(vcf_filepath):
    #generate_summary_sheet()
    #generate_gene_QC_sheet()
    #generate_exon_QC_sheet()
    generate_variant_sheets(vcf_filepath)


