import vcf2xls
import MySQLdb
import pytest
import vcf


class Test_geminiDB_query(object):
    def test_one(self):
        query = """ SELECT gid
                    FROM gene
                    WHERE name = 'CFTR'
                """
        
        assert vcf2xls.geminiDB_query(query) == [(5533L,)]


class Test_varfreq_query(object):
    def test_one(self):
        query = """ SELECT vid
                    FROM variant
                    WHERE chr = '1'
                    AND pos = '13372'
                    AND ref = 'G'
                    AND alt = 'C'
                """
        
        assert vcf2xls.varfreq_query(query) == [(1L,)]


class Test_execute_query(object):
    def test_one(self):
        query = """ SELECT gid
                    FROM gene
                    WHERE name = 'CFTR'
                """
        db = MySQLdb.connect(host="mgsrv01",
                         user="easih_ro",
                         db="GeminiDB")

        assert vcf2xls.execute_query(query, db) == [(5533L,)]


class Test_get_varfreq_vids(object):
    def test_one(self):
        chrom   = "1"
        pos     = "13372"
        ref     = "G"
        alt     = "C"

        assert vcf2xls.get_varfreq_vids(chrom, pos, ref, alt) == [1L]


class nonTest_gemini_allele_count(object):
    def test_one(self):
        # Can't be tested - output changes over time!
        '''
        with open("/data/projects/matt/vcf2xls/test/test_1_line.annotated.vcf") as vcf_fh:
            vcf_reader = vcf.Reader(vcf_fh)
        vcf_record = vcf_reader.next()
        alt        = G
        assert vcf2xls.gemini_allele_count(vcf_record, alt) ==\
            {"total":,
             "hom":,
             "het":
            }
        '''
        pass


class Test_sry_depth(object):
    def test_one(self):
        assert vcf2xls.sry_depth("G001010") == 114.25


class Test_get_read_stats(object):
    def test_one(self):
        assert vcf2xls.get_read_stats("G001010") == \
             {"total":31216116,
              "mapped":30925209,
              "duplicate":4116926}


class Test_gnum2pids(object):
    def test_one(self):
        g_number = "G002020"
        samplepanels_filepath = "/data/gemini/BioinformaticManifest.txt"
        assert vcf2xls.gnum2pids(g_number, samplepanels_filepath) == ["8969"]


class Test_pids2panelnames(object):
    def test_one(self):
        pids = ["8979","13"]
        genepanels_filepath = "/data/gemini/genepanels"
        assert vcf2xls.pids2panelnames(pids, genepanels_filepath) == \
            {"8979":"Achalasia-addisonianism-alacrimia syndrome",
             "13":"Xeroderma pigmentosum"}


class Test_pid2panelname(object):
    def test_one(self):
        pid = "8979"
        genepanels_filepath = "/data/gemini/genepanels"
        assert vcf2xls.pid2panelname(pid, genepanels_filepath) == "Achalasia-addisonianism-alacrimia syndrome"

    def test_two(self):
        pid = "_CFTR"
        genepanels_filepath = "/data/gemini/genepanels"
        assert vcf2xls.pid2panelname(pid, genepanels_filepath) == pid


class Test_pids2genes(object):
    def test_one(self):
        pids = ["8979","13"]
        genepanels_filepath = "/data/gemini/genepanels"
        assert vcf2xls.pids2genes(pids, genepanels_filepath) == \
            {"8979":["AAAS"],
             "13":  ["MPLKIP",
                     "GTF2H5",
                     "ERCC1",
                     "ERCC8",
                     "ERCC6",
                     "XPA",
                     "ERCC3",
                     "XPC",
                     "ERCC2",
                     "DDB2",
                     "ERCC4",
                     "ERCC5",
                     "POLH",
                    ]}


class Test_pid2genes(object):
    def test_one(self):
        pid = "13"
        genepanels_filepath = "/data/gemini/genepanels"
        assert vcf2xls.pid2genes(pid, genepanels_filepath) ==\
            ["MPLKIP",
             "GTF2H5",
             "ERCC1",
             "ERCC8",
             "ERCC6",
             "XPA",
             "ERCC3",
             "XPC",
             "ERCC2",
             "DDB2",
             "ERCC4",
             "ERCC5",
             "POLH",
                    ]


class Test_genes2transcripts(object):
    def test_one(self):
        genes = ["AAAS", "AARS", "AASS", "ARSE"]
        genes2transcripts_filepath = "/data/gemini/genes2transcripts"
        assert vcf2xls.genes2transcripts(genes, genes2transcripts_filepath) ==\
            {"AAAS":["NM_015665.5"],
             "AARS":["NM_001605.2"],
             "AASS":["NM_005763.3"],
             "ARSE":["NM_000047.2"],}


class Test_gene2transcript(object):
    def test_one(self):
        gene = "ARSE"
        genes2transcripts_filepath = "/data/gemini/genes2transcripts"
        assert vcf2xls.gene2transcripts(gene, genes2transcripts_filepath) ==\
            ["NM_000047.2"]


class Test_transcripts2exons(object):
    def test_one(self):
        transcripts = ["NM_001605.2","NM_000047.2"]
        assert vcf2xls.transcripts2exons(transcripts) ==\
            {"NM_001605.2":{"2":{   "chr":  "16",
                                    "start":70316523,
                                    "end":  70316666,},
                            "3":{   "chr":  "16",
                                    "start":70310869,
                                    "end":  70311057,},
                            "4":{   "chr":  "16",
                                    "start":70310389,
                                    "end":  70310534,},
                            "5":{   "chr":  "16",
                                    "start":70305684,
                                    "end":  70305875,},
                            "6":{   "chr":  "16",
                                    "start":70304099,
                                    "end":  70304243,},
                            "7":{   "chr":  "16",
                                    "start":70303521,
                                    "end":  70303666,},
                            "8":{   "chr":  "16",
                                    "start":70302174,
                                    "end":  70302282,},
                            "9":{   "chr":  "16",
                                    "start":70301562,
                                    "end":  70301712,},
                            "10":{  "chr":  "16",
                                    "start":70299441,
                                    "end":  70299565,},
                            "11":{  "chr":  "16",
                                    "start":70298861,
                                    "end":  70299005,},
                            "12":{  "chr":  "16",
                                    "start":70296249,
                                    "end":  70296427,},
                            "13":{  "chr":  "16",
                                    "start":70294947,
                                    "end":  70295060,},
                            "14":{  "chr":  "16",
                                    "start":70292883,
                                    "end":  70293089,},
                            "15":{  "chr":  "16",
                                    "start":70291936,
                                    "end":  70292120,},
                            "16":{  "chr":  "16",
                                    "start":70289631,
                                    "end":  70289739,},
                            "17":{  "chr":  "16",
                                    "start":70288524,
                                    "end":  70288637,},
                            "18":{  "chr":  "16",
                                    "start":70287822,
                                    "end":  70287941,},
                            "19":{  "chr":  "16",
                                    "start":70287617,
                                    "end":  70287703,},
                            "20":{  "chr":  "16",
                                    "start":70287171,
                                    "end":  70287284,},
                            "21":{  "chr":  "16",
                                    "start":70286624,
                                    "end":  70286809,},
                          },
            "NM_000047.2":{"2":{    "chr":  "X",
                                    "start":2878419,
                                    "end":  2878441,},
                            "3":{   "chr":  "X",
                                    "start":2876315,
                                    "end":  2876476,},
                            "4":{   "chr":  "X",
                                    "start":2873457,
                                    "end":  2873578,},
                            "5":{   "chr":  "X",
                                    "start":2871184,
                                    "end":  2871306,},
                            "6":{   "chr":  "X",
                                    "start":2867345,
                                    "end":  2867768,},
                            "7":{   "chr":  "X",
                                    "start":2864039,
                                    "end":  2864175,},
                            "8":{   "chr":  "X",
                                    "start":2861106,
                                    "end":  2861240,},
                            "9":{   "chr":  "X",
                                    "start":2856136,
                                    "end":  2856298,},
                            "10":{  "chr":  "X",
                                    "start":2854783,
                                    "end":  2854904,},
                            "11":{  "chr":  "X",
                                    "start":2852873,
                                    "end":  2853231,},
                        }
            }


class Test_transcripts2exons(object):
    def test_one(self):
        transcript = "NM_000047.2"
        assert vcf2xls.transcript2exons(transcript) ==\
                            {"2":{  "chr":  "X",
                                    "start":2878419,
                                    "end":  2878441,},
                            "3":{   "chr":  "X",
                                    "start":2876315,
                                    "end":  2876476,},
                            "4":{   "chr":  "X",
                                    "start":2873457,
                                    "end":  2873578,},
                            "5":{   "chr":  "X",
                                    "start":2871184,
                                    "end":  2871306,},
                            "6":{   "chr":  "X",
                                    "start":2867345,
                                    "end":  2867768,},
                            "7":{   "chr":  "X",
                                    "start":2864039,
                                    "end":  2864175,},
                            "8":{   "chr":  "X",
                                    "start":2861106,
                                    "end":  2861240,},
                            "9":{   "chr":  "X",
                                    "start":2856136,
                                    "end":  2856298,},
                            "10":{  "chr":  "X",
                                    "start":2854783,
                                    "end":  2854904,},
                            "11":{  "chr":  "X",
                                    "start":2852873,
                                    "end":  2853231,},
                        }


class Test_record_in_regions(object):
    def test_one(self):
        with open("/data/projects/matt/vcf2xls_py/test/test_1_line_missense.annotated.vcf") as vcf_fh:
            vcf_reader = vcf.Reader(vcf_fh)
            vcf_record = vcf_reader.next()
            gene_data_dict  = {"ARSE":{"transcripts":{"NM_000047.2":{"exons":{  "9":{   "chr":  "X",
                                                                                        "start": 2856136,
                                                                                        "end":   2856298,
                                                                                    },
                                                                                "10":{  "chr":  "X",
                                                                                        "start": 2854783,
                                                                                        "end":   2854904,
                                                                                    },

                                                                             }
                                                                    }
                                                    }
                                     }
                           }
            flank           = 30
        
            assert vcf2xls.record_in_regions(vcf_record, gene_data_dict, flank) == \
            ["NM_000047.2"]

    def test_two(self):
        with open("/data/projects/matt/vcf2xls_py/test/test_1_line_missense.annotated.vcf") as vcf_fh:
            vcf_reader = vcf.Reader(vcf_fh)
            vcf_record = vcf_reader.next()
            gene_data_dict  = {"ARSE":{"transcripts":{"NM_000047.2":{"exons":{  "11":{   "chr":  "X",
                                                                                        "start": 2852873,
                                                                                        "end":   2853231,
                                                                                    },
                                                                                "10":{  "chr":  "X",
                                                                                        "start": 2854783,
                                                                                        "end":   2854904,
                                                                                    },

                                                                             }
                                                                    }
                                                    }
                                     }
                           }
            flank           = 30
        
            assert vcf2xls.record_in_regions(vcf_record, gene_data_dict, flank) == \
            []    


class Test_record_range(object):
    def test_one(self):
        with open("/data/projects/matt/vcf2xls_py/test/test_1_line_missense.annotated.vcf") as vcf_fh:
            vcf_reader = vcf.Reader(vcf_fh)
            vcf_record = vcf_reader.next()
            assert vcf2xls.record_range(vcf_record) == [2856155, 2856155]


class Test_overlap(object):
    def test_one(self):
        # Full same multi
        range_one = [1, 10]
        range_two = [1, 10]
        assert vcf2xls.overlap(range_one, range_two) == 10

    def test_two(self):
        # One left
        range_one = [0, 10]
        range_two = [0, 0]
        assert vcf2xls.overlap(range_one, range_two) == 1

    def test_three(self):
        # Multi left
        range_one = [0, 10]
        range_two = [0, 1]
        assert vcf2xls.overlap(range_one, range_two) == 2

    def test_four(self):
        # One right
        range_one = [0, 10]
        range_two = [10, 10]
        assert vcf2xls.overlap(range_one, range_two) == 1

    def test_five(self):
        # Multi right
        range_one = [0, 10]
        range_two = [9, 10]
        assert vcf2xls.overlap(range_one, range_two) == 2

    def test_six(self):
        # One mid
        range_one = [0, 10]
        range_two = [5, 5]
        assert vcf2xls.overlap(range_one, range_two) == 1

    def test_seven(self):
        # Multi mid
        range_one = [0, 10]
        range_two = [5, 6]
        assert vcf2xls.overlap(range_one, range_two) == 2

    def test_eight(self):
        # Full same one
        range_one = [1, 1]
        range_two = [1, 1]
        assert vcf2xls.overlap(range_one, range_two) == 1

    def test_nine(self):
        # Adjacent none
        range_one = [1, 10]
        range_two = [11, 20]
        assert vcf2xls.overlap(range_one, range_two) == 0

    def test_ten(self):
        # Separate none
        range_one = [1, 10]
        range_two = [100, 110]
        assert vcf2xls.overlap(range_one, range_two) == 0

    def test_eleven(self):
        # Overlapping ends one
        range_one = [0, 10]
        range_two = [10, 20]
        assert vcf2xls.overlap(range_one, range_two) == 1

    def test_twelve(self):
        # Overlapping ends multi
        range_one = [0, 10]
        range_two = [9, 20]
        assert vcf2xls.overlap(range_one, range_two) == 2

    def test_thirteen(self):
        # No negatives
        with pytest.raises(AssertionError):
            range_one = [0, -1]
            range_two = [0, 20]
            vcf2xls.overlap(range_one, range_two)


class Test_get_csq_field_keys(object):
    def test_one(self):
        vcf_filepath = "/data/projects/matt/vcf2xls_py/test/test_1_line_missense.annotated.vcf"
        assert vcf2xls.get_csq_field_keys(vcf_filepath) == ["Allele","Gene","HGNC","RefSeq","Feature","Consequence","cDNA_position","Protein_position","Amino_acids","Existing_variation","SIFT","PolyPhen","HGVSc"]


class Test_chromosome_col(object):
    def test_one(self):
        alt = "Will be deleted"
        transcript = "Will be deleted"
        csq = "Will be deleted"
        csq_field_keys = "Will be deleted"
        with open("/data/projects/matt/vcf2xls_py/test/test_1_line_missense.annotated.vcf") as vcf_fh:
            vcf_reader = vcf.Reader(vcf_fh)
            vcf_record = vcf_reader.next()
            assert vcf2xls.chromosome_col(vcf_record, alt, transcript, csq, csq_field_keys) == \
                "X"


class nonTest_start_col(object):
    def test_one(self):
        alt = "Will be deleted"
        transcript = "Will be deleted"
        csq = "Will be deleted"
        csq_field_keys = "Will be deleted"
        with open("/data/projects/matt/vcf2xls_py/test/test_1_line_missense.annotated.vcf") as vcf_fh:
            vcf_reader = vcf.Reader(vcf_fh)
            vcf_record = vcf_reader.next()
            assert vcf2xls.start_col(vcf_record, alt, transcript, csq, csq_field_keys) == \
                "This is tricky - returns an excel formula object"


class Test_ref_col(object):
    def test_one(self):
        alt = "Will be deleted"
        transcript = "Will be deleted"
        csq = "Will be deleted"
        csq_field_keys = "Will be deleted"
        with open("/data/projects/matt/vcf2xls_py/test/test_1_line_missense.annotated.vcf") as vcf_fh:
            vcf_reader = vcf.Reader(vcf_fh)
            vcf_record = vcf_reader.next()
            assert vcf2xls.ref_col(vcf_record, alt, transcript, csq, csq_field_keys) == \
                "C"


class nonTest_pos_col(object):
    def test_one(self):
        alt = "Will be deleted"
        transcript = "Will be deleted"
        csq = "Will be deleted"
        csq_field_keys = "Will be deleted"
        with open("/data/projects/matt/vcf2xls_py/test/test_1_line_missense.annotated.vcf") as vcf_fh:
            vcf_reader = vcf.Reader(vcf_fh)
            vcf_record = vcf_reader.next()
            assert vcf2xls.pos_col(vcf_record, alt, transcript, csq, csq_field_keys) == \
                "This is tricky - returns an excel formula object"


class Test_qual_score_col(object):
    def test_one(self):
        alt = "Will be deleted"
        transcript = "Will be deleted"
        csq = "Will be deleted"
        csq_field_keys = "Will be deleted"
        with open("/data/projects/matt/vcf2xls_py/test/test_1_line_missense.annotated.vcf") as vcf_fh:
            vcf_reader = vcf.Reader(vcf_fh)
            vcf_record = vcf_reader.next()
            assert vcf2xls.qual_score_col(vcf_record, alt, transcript, csq, csq_field_keys) == \
                1868.77


class Test_aaf_col(object):
    def test_one(self):
        alt = "Will be deleted"
        transcript = "Will be deleted"
        csq = "Will be deleted"
        csq_field_keys = "Will be deleted"
        with open("/data/projects/matt/vcf2xls_py/test/test_1_line_missense.annotated.vcf") as vcf_fh:
            vcf_reader = vcf.Reader(vcf_fh)
            vcf_record = vcf_reader.next()
            assert vcf2xls.aaf_col(vcf_record, alt, transcript, csq, csq_field_keys) ==\
                0.429


class Test_genotype_col(object):
    def test_one(self):
        alt = "Will be deleted"
        transcript = "Will be deleted"
        csq = "Will be deleted"
        csq_field_keys = "Will be deleted"
        with open("/data/projects/matt/vcf2xls_py/test/test_1_line_missense.annotated.vcf") as vcf_fh:
            vcf_reader = vcf.Reader(vcf_fh)
            vcf_record = vcf_reader.next()
            assert vcf2xls.genotype_col(vcf_record, alt, transcript, csq, csq_field_keys) ==\
                ["0/1"]


class nonTest_dbsnp_col(object):
    def test_one(self):
        alt = "Will be deleted"
        transcript = "Will be deleted"
        csq = "Will be deleted"
        csq_field_keys = "Will be deleted"
        with open("/data/projects/matt/vcf2xls_py/test/test_1_line_missense.annotated.vcf") as vcf_fh:
            vcf_reader = vcf.Reader(vcf_fh)
            vcf_record = vcf_reader.next()
            assert vcf2xls.dbsnp_col(vcf_record, alt, transcript, csq, csq_field_keys) ==\
                "This is tricky - returns an excel formula object"


class nonTest_af_gemini_col(object):
    def test_one(self):
        alt = "Will be deleted"
        transcript = "Will be deleted"
        csq = "Will be deleted"
        csq_field_keys = "Will be deleted"
        with open("/data/projects/matt/vcf2xls_py/test/test_1_line_missense.annotated.vcf") as vcf_fh:
            vcf_reader = vcf.Reader(vcf_fh)
            vcf_record = vcf_reader.next()
            assert vcf2xls.af_gemini_col(vcf_record, alt, transcript, csq, csq_field_keys) ==\
                "This is tricky - value changes over time"


class nonTest_het_hom_gemini_col(object):
    def test_one(self):
        alt = "Will be deleted"
        transcript = "Will be deleted"
        csq = "Will be deleted"
        csq_field_keys = "Will be deleted"
        with open("/data/projects/matt/vcf2xls_py/test/test_1_line_missense.annotated.vcf") as vcf_fh:
            vcf_reader = vcf.Reader(vcf_fh)
            vcf_record = vcf_reader.next()
            assert vcf2xls.het_hom_gemini_col(vcf_record, alt, transcript, csq, csq_field_keys) ==\
                "This is tricky - value changes over time"


class nonTest_exac_col(object):
    def test_one(self):
        alt = "T"
        transcript = "Will be deleted"
        csq = "Will be deleted"
        csq_field_keys = "Will be deleted"
        with open("/data/projects/matt/vcf2xls_py/test/test_1_line_missense.annotated.vcf") as vcf_fh:
            vcf_reader = vcf.Reader(vcf_fh)
            vcf_record = vcf_reader.next()
            assert vcf2xls.af_exac_col(vcf_record, alt, transcript, csq, csq_field_keys) ==\
                "0.6593" # http://exac.broadinstitute.org/variant/X-2856155-C-T

class Test_depth_col(object):
    def test_one(self):
        alt = "Will be deleted"
        transcript = "Will be deleted"
        csq = "Will be deleted"
        csq_field_keys = "Will be deleted"
        with open("/data/projects/matt/vcf2xls_py/test/test_1_line_missense.annotated.vcf") as vcf_fh:
            vcf_reader = vcf.Reader(vcf_fh)
            vcf_record = vcf_reader.next()
            assert vcf2xls.depth_col(vcf_record, alt, transcript, csq, csq_field_keys) ==\
                ["133"]

class Test_gene_col(object):
    def test_one(self):
        alt = "Will be deleted"
        transcript = "Will be deleted"
        csq = "T|ENSG00000157399|ARSE|NM_000047.2|ENST00000381134|missense_variant|1337|424|G/S||tolerated(0.21)|benign(0.132)|ENST00000381134.3:c.1270G>A"
        csq_field_keys = ["Allele","Gene","HGNC","RefSeq","Feature","Consequence","cDNA_position","Protein_position","Amino_acids","Existing_variation","SIFT","PolyPhen","HGVSc"]
        vcf_record = "Will be deleted"
        assert vcf2xls.gene_col(vcf_record, alt, transcript, csq, csq_field_keys) ==\
            "ARSE"

class Test_c_pos_col(object):
    def test_one(self):
        alt = "Will be deleted"
        transcript = "Will be deleted"
        csq = "T|ENSG00000157399|ARSE|NM_000047.2|ENST00000381134|missense_variant|1337|424|G/S||tolerated(0.21)|benign(0.132)|ENST00000381134.3:c.1270G>A"
        csq_field_keys = ["Allele","Gene","HGNC","RefSeq","Feature","Consequence","cDNA_position","Protein_position","Amino_acids","Existing_variation","SIFT","PolyPhen","HGVSc"]
        vcf_record = "Will be deleted"
        assert vcf2xls.c_pos_col(vcf_record, alt, transcript, csq, csq_field_keys) ==\
            "c.1270"


class Test_c_change_col(object):
    def test_one(self):
        alt = "Will be deleted"
        transcript = "Will be deleted"
        csq = "T|ENSG00000157399|ARSE|NM_000047.2|ENST00000381134|missense_variant|1337|424|G/S||tolerated(0.21)|benign(0.132)|ENST00000381134.3:c.1270G>A"
        csq_field_keys = ["Allele","Gene","HGNC","RefSeq","Feature","Consequence","cDNA_position","Protein_position","Amino_acids","Existing_variation","SIFT","PolyPhen","HGVSc"]
        vcf_record = "Will be deleted"
        assert vcf2xls.c_change_col(vcf_record, alt, transcript, csq, csq_field_keys) ==\
            "G>A"

class Test_aa_change_col(object):
    def test_one(self):
        alt = "Will be deleted"
        transcript = "Will be deleted"
        csq = "T|ENSG00000157399|ARSE|NM_000047.2|ENST00000381134|missense_variant|1337|424|G/S||tolerated(0.21)|benign(0.132)|ENST00000381134.3:c.1270G>A"
        csq_field_keys = ["Allele","Gene","HGNC","RefSeq","Feature","Consequence","cDNA_position","Protein_position","Amino_acids","Existing_variation","SIFT","PolyPhen","HGVSc"]
        vcf_record = "Will be deleted"
        assert vcf2xls.aa_change_col(vcf_record, alt, transcript, csq, csq_field_keys) ==\
            "G424S"

class Test_consequences_col(object):
    def test_one(self):
        alt = "Will be deleted"
        transcript = "Will be deleted"
        csq = "T|ENSG00000157399|ARSE|NM_000047.2|ENST00000381134|missense_variant|1337|424|G/S||tolerated(0.21)|benign(0.132)|ENST00000381134.3:c.1270G>A"
        csq_field_keys = ["Allele","Gene","HGNC","RefSeq","Feature","Consequence","cDNA_position","Protein_position","Amino_acids","Existing_variation","SIFT","PolyPhen","HGVSc"]
        vcf_record = "Will be deleted"
        assert vcf2xls.consequences_col(vcf_record, alt, transcript, csq, csq_field_keys) ==\
            "missense_variant"

class Test_sift_col(object):
    def test_one(self):
        alt = "Will be deleted"
        transcript = "Will be deleted"
        csq = "T|ENSG00000157399|ARSE|NM_000047.2|ENST00000381134|missense_variant|1337|424|G/S||tolerated(0.21)|benign(0.132)|ENST00000381134.3:c.1270G>A"
        csq_field_keys = ["Allele","Gene","HGNC","RefSeq","Feature","Consequence","cDNA_position","Protein_position","Amino_acids","Existing_variation","SIFT","PolyPhen","HGVSc"]
        vcf_record = "Will be deleted"
        assert vcf2xls.sift_col(vcf_record, alt, transcript, csq, csq_field_keys) ==\
            "tolerated(0.21)"

class Test_polyphen_col(object):
    def test_one(self):
        alt = "Will be deleted"
        transcript = "Will be deleted"
        csq = "T|ENSG00000157399|ARSE|NM_000047.2|ENST00000381134|missense_variant|1337|424|G/S||tolerated(0.21)|benign(0.132)|ENST00000381134.3:c.1270G>A"
        csq_field_keys = ["Allele","Gene","HGNC","RefSeq","Feature","Consequence","cDNA_position","Protein_position","Amino_acids","Existing_variation","SIFT","PolyPhen","HGVSc"]
        vcf_record = "Will be deleted"
        assert vcf2xls.polyphen_col(vcf_record, alt, transcript, csq, csq_field_keys) ==\
            "benign(0.132)" 


class nonTest_OMIM_col(object):
    def test_one(self):
        alt = "Will be deleted"
        transcript = "Will be deleted"
        csq = "T|ENSG00000157399|ARSE|NM_000047.2|ENST00000381134|missense_variant|1337|424|G/S||tolerated(0.21)|benign(0.132)|ENST00000381134.3:c.1270G>A"
        csq_field_keys = ["Allele","Gene","HGNC","RefSeq","Feature","Consequence","cDNA_position","Protein_position","Amino_acids","Existing_variation","SIFT","PolyPhen","HGVSc"]
        vcf_record = "Will be deleted"
        assert vcf2xls.OMIM_col(vcf_record, alt, transcript, csq, csq_field_keys) ==\
            "Tricky"


class Test_transcript_col(object):
    def test_one(self):
        alt = "Will be deleted"
        transcript = "NM_000047.2"
        csq = "T|ENSG00000157399|ARSE|NM_000047.2|ENST00000381134|missense_variant|1337|424|G/S||tolerated(0.21)|benign(0.132)|ENST00000381134.3:c.1270G>A"
        csq_field_keys = ["Allele","Gene","HGNC","RefSeq","Feature","Consequence","cDNA_position","Protein_position","Amino_acids","Existing_variation","SIFT","PolyPhen","HGVSc"]
        vcf_record = "Will be deleted"
        assert vcf2xls.transcript_col(vcf_record, alt, transcript, csq, csq_field_keys) ==\
            "NM_000047.2"


class nonTest_alt_col(object):
    def test_one(self):
        alt = "Will be deleted"
        transcript = "Will be deleted"
        csq = "T|ENSG00000157399|ARSE|NM_000047.2|ENST00000381134|missense_variant|1337|424|G/S||tolerated(0.21)|benign(0.132)|ENST00000381134.3:c.1270G>A"
        csq_field_keys = ["Allele","Gene","HGNC","RefSeq","Feature","Consequence","cDNA_position","Protein_position","Amino_acids","Existing_variation","SIFT","PolyPhen","HGVSc"]
        vcf_record = "Will be deleted"
        assert vcf2xls.alt_col(vcf_record, alt, transcript, csq, csq_field_keys) ==\
            "Trivial"

class nonTest_af_gnomad_exome_col(object):
    def test_one(self):
        alt = "Will be deleted"
        transcript = "Will be deleted"
        csq = "T|ENSG00000157399|ARSE|NM_000047.2|ENST00000381134|missense_variant|1337|424|G/S||tolerated(0.21)|benign(0.132)|ENST00000381134.3:c.1270G>A"
        csq_field_keys = ["Allele","Gene","HGNC","RefSeq","Feature","Consequence","cDNA_position","Protein_position","Amino_acids","Existing_variation","SIFT","PolyPhen","HGVSc"]
        vcf_record = "Will be deleted"
        assert vcf2xls.af_gnomad_exome_col(vcf_record, alt, transcript, csq, csq_field_keys) ==\
            "Tricky"

class nonTest_af_gnomad_genome_col(object):
    def test_one(self):
        alt = "Will be deleted"
        transcript = "Will be deleted"
        csq = "T|ENSG00000157399|ARSE|NM_000047.2|ENST00000381134|missense_variant|1337|424|G/S||tolerated(0.21)|benign(0.132)|ENST00000381134.3:c.1270G>A"
        csq_field_keys = ["Allele","Gene","HGNC","RefSeq","Feature","Consequence","cDNA_position","Protein_position","Amino_acids","Existing_variation","SIFT","PolyPhen","HGVSc"]
        vcf_record = "Will be deleted"
        assert vcf2xls.af_gnomad_genome_col(vcf_record, alt, transcript, csq, csq_field_keys) ==\
            "Tricky"

class Test_get_csq_field(object):
    def test_one(self):
        field_name = "Feature"
        csq = "T|ENSG00000157399|ARSE|NM_000047.2|ENST00000381134|missense_variant|1337|424|G/S||tolerated(0.21)|benign(0.132)|ENST00000381134.3:c.1270G>A"
        csq_field_keys = ["Allele","Gene","HGNC","RefSeq","Feature","Consequence","cDNA_position","Protein_position","Amino_acids","Existing_variation","SIFT","PolyPhen","HGVSc"]
        assert vcf2xls.get_csq_field(field_name, csq, csq_field_keys) ==\
            "ENST00000381134"

class Test_alts2csqallele(object):
    def test_one(self):
        with open("/data/projects/matt/vcf2xls_py/test/test_1_line_missense.annotated.vcf") as vcf_fh:
            vcf_reader = vcf.Reader(vcf_fh)
            vcf_record = vcf_reader.next()
        assert vcf2xls.alts2csqallele(vcf_record) ==\
            ["T"]





