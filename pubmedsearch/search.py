# incorporates search functions for medline and entrez using biopython

from Bio import Entrez, Medline


class Search:

    def __init__(self, email, query, idlist):
        self.email = email
        self.query = query
        self.idlist = idlist
        self.count = []

    def search_entrez(self):
        Entrez.email = self.email
        handle = Entrez.egquery(term=self.query)
        record_0 = Entrez.read(handle)
        for row in record_0["eGQueryResult"]:
            if row["DbName"] == "pubmed":
                self.count = (row["Count"])

        print("\n", "Your search query returned", self.count, "results from PubMed.")
        print("\n", "Now retrieving records from MedLine...")

        handle_1 = Entrez.esearch(db="pubmed", term=self.query, retmax=self.count)
        record_1 = Entrez.read(handle_1)
        handle_1.close()
        self.idlist = record_1["IdList"]

        return self.idlist

    def search_medline(self):
        handle_2 = Entrez.efetch(db="pubmed", id=self.idlist, rettype="medline", retmode="text")
        record_2 = Medline.parse(handle_2)
        records = list(record_2)
        print("Done")
        print("\n", len(records), "records returned from Medline.")

        return records
