import json
from collections import OrderedDict
from multiqc.modules.base_module import BaseMultiqcModule
from multiqc import config
from multiqc.plots import bar

class MultiqcModule(BaseMultiqcModule):
    def __init__(self):
        # 1. Initialise the parent class
        super(MultiqcModule, self).__init__(
            name='RNA-Seq',
            anchor='multiomic_report',
            href='https://github.com/regulatory-genomics/RNA-sm',
            info='custom pipeline'
        )

        # -----------------------------------------------------------
        # 2. DATA STORAGE
        # -----------------------------------------------------------
        self.rnaseqc_data = dict()    # Stores metrics.tsv data
        self.genetype_data = dict()   # Stores gene_type_counts data
        self.rsem_data = dict()       # Stores rsem.genes.json data

        # -----------------------------------------------------------
        # 3. PARSING LOGIC
        # -----------------------------------------------------------
        
        # A. Parse RNA-SeQC (TSV)
        # Matches 'rnaseq/rnaseqqc' from your execution_start
        for f in self.find_log_files('rnaseq/rnaseqqc'):
            self.parse_rnaseqc_tsv(f)

        # B. Parse Gene Types (JSON)
        # Matches 'rnaseq/gene_type_counts' from your execution_start
        for f in self.find_log_files('rnaseq/gene_type_counts'):
            self.parse_genetype_json(f)

        # C. Parse RSEM Genes (JSON)
        # Matches 'rnaseq/rsem' from your execution_start
        for f in self.find_log_files('rnaseq/rsem'):
            self.parse_rsem_json(f)

        # -----------------------------------------------------------
        # 4. FILTERING & EXIT
        # -----------------------------------------------------------
        # Filter out samples ignored by the user
        self.rnaseqc_data = self.ignore_samples(self.rnaseqc_data)
        self.genetype_data = self.ignore_samples(self.genetype_data)
        self.rsem_data = self.ignore_samples(self.rsem_data)

        # If no data found at all, raise warning
        if len(self.rnaseqc_data) == 0 and len(self.genetype_data) == 0:
            raise UserWarning

        # -----------------------------------------------------------
        # 5. GENERATE REPORT SECTIONS
        # -----------------------------------------------------------
        self.write_general_stats()
        self.write_gene_type_plot()

    # ===============================================================
    # PARSING FUNCTIONS
    # ===============================================================

    def parse_rnaseqc_tsv(self, f):
        """ Parses the RNA-SeQC TSV file (Key [tab] Value) """
        parsed_data = {}
        for line in f['content'].splitlines():
            s = line.split('\t')
            if len(s) > 1:
                key = s[0].strip()
                val = s[1].strip()
                # Try to convert numbers to floats
                try:
                    parsed_data[key] = float(val)
                except ValueError:
                    parsed_data[key] = val # Keep strings (like Sample name)
        
        if parsed_data:
            self.rnaseqc_data[f['s_name']] = parsed_data

    def parse_genetype_json(self, f):
        """ Parses the nested Gene Type JSON file """
        try:
            data = json.loads(f['content'])
            # The file structure is {"gene_type_count": {...}}
            if "gene_type_count" in data:
                self.genetype_data[f['s_name']] = data["gene_type_count"]
        except Exception as e:
            pass # Skip malformed files

    def parse_rsem_json(self, f):
        """ Parses the RSEM JSON file """
        try:
            data = json.loads(f['content'])
            if "num_genes_detected" in data:
                self.rsem_data[f['s_name']] = data
        except Exception as e:
            pass

    # ===============================================================
    # REPORT WRITING
    # ===============================================================

    def write_general_stats(self):
        """ Adds columns to the main top table in MultiQC """
        
        # 1. Define Headers for RNA-SeQC Metrics
        rnaseqc_headers = OrderedDict()
        rnaseqc_headers['Mapping Rate'] = {
            'title': 'Map Rate',
            'description': 'RNA-SeQC: Mapping Rate',
            'max': 1, 'min': 0, 'suffix': '%',
            'scale': 'PuBu',
            'format': '{:,.1%}' # 0.41 -> 41.0%
        }
        rnaseqc_headers['Exonic Rate'] = {
            'title': 'Exonic',
            'description': 'RNA-SeQC: Exonic Rate',
            'max': 1, 'min': 0, 'suffix': '%',
            'scale': 'Greens',
            'format': '{:,.1%}'
        }
        
        # 2. Define Headers for RSEM Metrics
        rsem_headers = OrderedDict()
        rsem_headers['num_genes_detected'] = {
            'title': 'Genes',
            'description': 'RSEM: Number of Genes Detected',
            'format': '{:,.0f}',
            'scale': 'OrRd'
        }

        # 3. Add to General Stats
        # We call this twice to merge data from different dictionaries
        self.general_stats_addcols(self.rnaseqc_data, rnaseqc_headers)
        self.general_stats_addcols(self.rsem_data, rsem_headers)

    def write_gene_type_plot(self):
        """ Creates a Stacked Bar Plot for Gene Types """
        if not self.genetype_data:
            return

        # Configuration for the plot
        pconfig = {
            'id': 'my_rnaseq_genetypes',
            'title': 'My RNA-Seq: Gene Type Counts',
            'ylab': 'Count',
            'cpswitch_counts_label': 'Number of Genes'
        }
        
        # Sort categories so the legend is clean
        # (Optional: MultiQC sorts automatically, but this helps consistency)
        cats = sorted(list(set(k for s in self.genetype_data.values() for k in s)))
        
        self.add_section(
            name='Gene Types',
            anchor='my_rnaseq_genetypes',
            description='Counts of different gene biotypes (protein_coding, rRNA, etc.)',
            plot=bar.plot(self.genetype_data, cats, pconfig)
        )