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
    variant_types = ["de_novo_AD", "other"]
    
    # might need to make each a dict of key = name, value = list of seq ontology terms for type
    for variant_type in variant_types[0:1]:
        ws = wb.add_sheet(variant_type)
        generate_variant_sheet(variant_type, vcf_filepath, ws)


def generate_variant_sheet(variant_type, vcf_filepath, worksheet):
    vcf_reader = vcf.Reader(open(vcf_filepath, 'r'))
    csq_fields = get_CSQ_format(vcf_reader)
    sheet_columns = []
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

def is_AR(vcf_record):
    # Criteria:
    # Variant call present in proband
    # Proband genotype is hom alt (2 identical non-ref alleles)
    # Neither parent genotype is hom alt

    proband_call, maternal_call, paternal_call = vcf_record.samples
    
    # Some sites will only be called in parents
    if proband_call.called:
    
        # Proband is homozygous and genotype contains an alt allele
        # Therefore must be hom-alt
        if is_homozygous(proband_call) and contains_alt(proband_call):
            
            # Variant not hom-alt in either parent
            for parental_call in [maternal_call, paternal_call]:
                if is_same_genotype(proband_call.gt_alleles, parental_call):
                    return False
            return True
    return False


def is_AD(vcf_record):
    # Criteria:
    # Variant call present in proband
    # Proband genotype is het ref/alt
    # Neither parent genotype contains the alt

    proband_call, maternal_call, paternal_call = vcf_record.samples
    
    # Some sites will only be called in parents
    if proband_call.called:
    
        # Proband is heterozygous and genotype contains ref allele
        # Therefore must be het ref-alt (but not 1/2)
        if is_heterozygous(proband_call) and contains_ref(proband_call):
            
            # Variant should not be present in either parent
            # First find the ID of the alt allele in the proband (1/2/3 etc)
            proband_alt_alleles = get_alt_alleles(proband_call)
            assert len(proband_alt_alleles) == 1, "Het ref alt GT should contain only one alt"
            proband_alt_allele = proband_alt_alleles[0]
            
            # Search parental genotypes for that alt
            for parental_call in [maternal_call, paternal_call]:
                if contains_allele(proband_alt_allele, parental_call):
                    return False
            return True
    return False


def is_homozygous(call):
    # Criteria:
    # Both alleles in a call are the same

    assert len(call.gt_alleles) == 2, "Record does not contain 2 alleles"
    if call.gt_alleles[0] == call.gt_alleles[1]:
        return True
    else:
        return False


def is_heterozygous(call):
    # Criteria:
    # Both alleles in a call are different

    assert len(call.gt_alleles) == 2, "Record does not contain 2 alleles"
    if call.gt_alleles[0] != call.gt_alleles[1]:
        return True
    else:
        return False


def is_same_genotype(query_genotype, subject_call):
    # If subject allele is missing we cannot tell if they are the same
    if not subject_call.called:
        return False
    else:
        if query_genotype == subject_call.gt_alleles:
            
            return True
        else:
            return False


def contains_alt(call):
    # Criteria:
    # One or more alleles in a call non-ref

    alt_allele = any([True for x in call.gt_alleles if x != "0"])
    if alt_allele:
        return True
    else:
        return False


def contains_allele(query_allele, call):
    # Criteria:
    # Genotype of call contains query allele

    if call.called:
        # Any allele in call matches the query allele
        allele_found = query_allele in call.gt_alleles
        if allele_found:
            return True
        else:
            return False
    return False


def contains_ref(call):
    # Criteria:
    # One or more alleles in a call are ref

    ref_allele = any([True for x in call.gt_alleles if x == "0"])
    if ref_allele:
        return True
    else:
        return False


def get_alt_alleles(call):
    # Return a list of all alt alleles in the genotype
    alt_alleles = [x for x in call.gt_alleles if x != "0"]
    return alt_alleles


def is_denovo(vcf_record):

    family = {"G003906": "proband",
              "G003907": "mother",
              "G003908": "father",  }
    
    proband_call, maternal_call, paternal_call = vcf_record.samples
    
    if proband_call.called and not proband_call.phased:

        # At this point we have all unphased proband calls
        # so triple hets and de novos
        # also called proband + missing parent(s)

        parental_calls = \
            [call for call in [maternal_call, paternal_call] if call.called]
        
        # Non ref alleles only
        proband_variant_alleles = [x for x in proband_call.gt_alleles if x != 0]

        for proband_allele in proband_variant_alleles:
            # Start from the assumption that all calls are de-novo
            denovo = True
            
            # Look for conflicting or missing evidence
            # Missing parental calls - could be de novo
            if len(parental_calls) < 2:
                print vcf_record
                print vcf_record.samples
                print "MAYBE DE NOVO"
                raw_input()
                break

            else:
                for parental_call in parental_calls:
                    #print parental_call.gt_alleles
                    if proband_allele in parental_call.gt_alleles:
                        #print "NOT DE NOVO"
                        denovo = False
                        break
                
                if denovo == True:
                    print vcf_record
                    print vcf_record.samples
                    print "DE NOVO"
                    raw_input()

        #    parental_calls = [call for call in parental_calls if parental_calls.called]
        #    print parental_calls
        #    if proband_allele not in parental_alleles:
        #        print "DE NOVO"
        #if proband_maternal_allele not in maternal_call.gt_alleles:
        #            print "DE NOVO!"



    '''
    print family[call.sample]

    try:
        print call.gt_alleles
        print call.phased


        except AttributeError:
            print None
    '''


def get_next_region(regions_fh):
    line = regions_fh.readline()
    if line:
        chrom, start, end, name = line.strip().split()
        region = { "chrom": chrom,
                   "start": start,
                     "end": end,
                    "name": name
                 }
        return region
    return None

def coordinate_overlap(a, b):
    return max(0, min(a[1], b[1]) - max(a[0], b[0]))

def assign_variants(trio_phased_vcf):
    regions_fh = open("/data/projects/matt/trios/data/trio.bed")
    region = get_next_region(regions_fh)
    records_in_region = []

    vcf_fh = open_fh(trio_phased_vcf, read_write="r", ext=".vcf")
    vcf_reader = vcf.Reader(vcf_fh)
    

    while vcf_reader:
        try:
            current_row = vcf_reader.next()
        except StopIteration:
            print "\nEOF"
            break

        if current_row:
            # is the variant on the right chrom
            while current_row.CHROM > region["chrom"]:
                region = get_next_region(regions_fh)


            end = current_row.POS + max([len(x) for x in [current_row.REF] + current_row.ALT])
            print "\nVAR:"
            print current_row.CHROM, current_row.POS, end, current_row.REF, current_row.ALT
            print "Region:"
            print region["chrom"], region["start"], region["end"]
            if current_row.POS >= region["start"] and current_row.POS <= region["end"]
        else:
            break
            

    exit()

    #for row_index, vcf_record in enumerate(vcf_reader):
        # get rows until current_row falls outside current region
        # get regions until region pos = or > current_row

        # To do - store rows by region here
    '''
        if pos in region:
            add record to region list
        else:
            get next region until region pos same or greater than pos
            reset list





        for each region
            list of vars = []
            while vcf_reader contains stuff
                get lines and add them to the list until a line outside the list if found

    '''
    '''    
        # Autosomal recessive
        if is_AR(vcf_record):
            print
            print "AR"
            print vcf_record
            print vcf_record.samples

        # Autosomal dominant
        elif is_AD(vcf_record):
            print
            print "AD"
            print vcf_record
            print vcf_record.samples
        
        # The bin
        else:
            print
            print "Other"
            print vcf_record
            print vcf_record.samples
    '''

@begin.start
def main(vcf_filepath):
    #generate_summary_sheet()
    #generate_gene_QC_sheet()
    #generate_exon_QC_sheet()
    #generate_variant_sheets(vcf_filepath)
    assign_variants(vcf_filepath)

