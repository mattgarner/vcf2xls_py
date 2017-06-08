import vcf2xls
import MySQLdb


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

        assert get_varfreq_vids(chrom, pos, ref, alt) == [(1L,)]


class Test_gemini_allele_count(object):
    def test_one(self):
        vcf_record = 1   
        alt        = 1
        assert 1==1

