import pysam

vcf_filepath = "/data/gemini/170726_170728_NS500192/vcfs/X004877.annotated.vcf"
vcf_in = pysam.VariantFile(vcf_filepath)  # auto-detect input format

print vcf_in.__doc__
print vcf_in.header




'''
with open(vcf_filepath) as vcf_fh:
    vcf_reader = vcf.Reader(vcf_fh)
    csq_desc = vcf_reader.infos.get("CSQ").desc
    csq_field_keys = csq_desc.split("Format: ")[-1].split("|")
    return csq_field_keys
'''