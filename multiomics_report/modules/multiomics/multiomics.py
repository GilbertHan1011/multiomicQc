import json
import logging
from collections import OrderedDict
from multiqc.base_module import BaseMultiqcModule, ModuleNoSamplesFound
from multiqc import config
from multiqc.plots import bargraph

log = logging.getLogger(__name__)

class MultiqcModule(BaseMultiqcModule):
    def __init__(self):
        # 1. Initialise the parent class
        super(MultiqcModule, self).__init__(
            name='multiomics',
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
        self.mad_data = dict()        # Stores mad_qc.json data
        # -----------------------------------------------------------
        # 3. PARSING LOGIC
        # -----------------------------------------------------------
        
        # A. Parse RNA-SeQC (TSV)
        # Matches 'multiomics/rnaseqqc' from your execution_start
        for f in self.find_log_files('multiomics_report/rnaseqqc'):
            self.parse_rnaseqc_tsv(f)

        # B. Parse Gene Types (JSON)
        # Matches 'multiomics/gene_type_counts' from your execution_start
        for f in self.find_log_files('multiomics_report/gene_type_counts'):
            self.parse_genetype_json(f)

        # C. Parse RSEM Genes (JSON)
        # Matches 'multiomics/rsem' from your execution_start
        for f in self.find_log_files('multiomics_report/rsem'):
            self.parse_rsem_json(f)

        # Debug: Log search for MAD files
        log.info("[MAD DEBUG] Searching for files matching 'multiomics_report/mad_rna' pattern...")
        mad_files = list(self.find_log_files('multiomics_report/mad_rna'))
        log.info(f"[MAD DEBUG] Found {len(mad_files)} file(s) matching 'multiomics_report/mad_rna' pattern")
        
        if len(mad_files) == 0:
            # Additional debugging: check what files are available
            log.warning("[MAD DEBUG] No files found! This could mean:")
            log.warning("[MAD DEBUG]   1. No files match the regex pattern: '.*summary_metrics\\.json$'")
            log.warning("[MAD DEBUG]   2. Files match filename but don't contain 'MAD of log ratios' in first 10 lines")
            log.warning("[MAD DEBUG]   3. Files are in directories not searched by MultiQC")
        
        for f in mad_files:
            log.info(f"[MAD DEBUG] Processing MAD file: {f.get('fn', 'unknown')} (sample: {f.get('s_name', 'unknown')})")
            self.parse_mad_json(f)
        # -----------------------------------------------------------
        # 4. FILTERING & EXIT
        # -----------------------------------------------------------
        # Filter out samples ignored by the user
        self.rnaseqc_data = self.ignore_samples(self.rnaseqc_data)
        self.genetype_data = self.ignore_samples(self.genetype_data)
        self.rsem_data = self.ignore_samples(self.rsem_data)
        self.mad_data = self.ignore_samples(self.mad_data)
        # If no data found at all, raise ModuleNoSamplesFound
        if (len(self.rnaseqc_data) == 0 and len(self.genetype_data) == 0 and 
            len(self.rsem_data) == 0 and len(self.mad_data) == 0):
            raise ModuleNoSamplesFound

        # -----------------------------------------------------------
        # 5. GENERATE REPORT SECTIONS
        # -----------------------------------------------------------
        # Add software version (required by MultiQC)
        self.add_software_version(None)
        
        self.write_general_stats()
        self.write_gene_type_plot()
        
        # Write data file (MUST be at the end, after all sections are added)
        # Combine all data into one dict for writing
        all_data = {
            'rnaseqc': self.rnaseqc_data,
            'genetype': self.genetype_data,
            'rsem': self.rsem_data,
            'mad': self.mad_data
        }
        self.write_data_file(all_data, "multiqc_multiomics")

    # ===============================================================
    # PARSING FUNCTIONS
    # ===============================================================

    def parse_rnaseqc_tsv(self, f):
        """ Parses the RNA-SeQC TSV file (Key [tab] Value) """
        parsed_data = {}
        for line in f['f'].splitlines():
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
            data = json.loads(f['f'])
            # The file structure is {"gene_type_count": {...}}
            if "gene_type_count" in data:
                self.genetype_data[f['s_name']] = data["gene_type_count"]
        except Exception as e:
            pass # Skip malformed files

    def parse_rsem_json(self, f):
        """ Parses the RSEM JSON file """
        try:
            data = json.loads(f['f'])
            if "num_genes_detected" in data:
                self.rsem_data[f['s_name']] = data
        except Exception as e:
            pass
    
    def parse_mad_json(self, f):
        """ Parses the MAD QC summary JSON file """
        try:
            log.debug(f"[MAD DEBUG] Parsing file: {f.get('fn', 'unknown')}")
            file_content = f.get('f', '')
            log.debug(f"[MAD DEBUG] File content length: {len(file_content)} chars")
            log.debug(f"[MAD DEBUG] First 200 chars of content: {file_content[:200]}")
            
            data = json.loads(file_content)
            log.debug(f"[MAD DEBUG] JSON parsed successfully. Keys: {list(data.keys())[:10]}")
            
            # Validation: specific key check
            if "MAD of log ratios" in data:
                log.info(f"[MAD DEBUG] Found 'MAD of log ratios' key in file {f.get('fn', 'unknown')}")
                # Use the clean sample name provided by MultiQC
                s_name = f['s_name']
                
                # Check if this sample already has data (handling duplicates)
                if s_name in self.mad_data:
                    log.debug(f"[MAD DEBUG] Duplicate sample found for MAD QC: {s_name}")
                
                # Also check for status field to skip invalid files
                if "status" in data and data["status"] == "No pairs available for aggregation":
                    # Skip files with no pairs
                    log.warning(f"[MAD DEBUG] Skipping file {f.get('fn', 'unknown')}: No pairs available")
                    return
                
                self.mad_data[s_name] = data
                self.add_data_source(f, section='mad_qc')
                log.info(f"[MAD DEBUG] Successfully added MAD data for sample: {s_name}")
            else:
                log.warning(f"[MAD DEBUG] File {f.get('fn', 'unknown')} does not contain 'MAD of log ratios' key")
                log.debug(f"[MAD DEBUG] Available keys in JSON: {list(data.keys())}")
                
        except json.JSONDecodeError as e:
            log.warning(f"[MAD DEBUG] Failed to parse MAD QC JSON file {f.get('fn', 'unknown')}: {e}")
            log.debug(f"[MAD DEBUG] Content preview: {f.get('f', '')[:500]}")
        except Exception as e:
            log.warning(f"[MAD DEBUG] Error parsing MAD QC file {f.get('fn', 'unknown')}: {e}")
            import traceback
            log.debug(f"[MAD DEBUG] Traceback: {traceback.format_exc()}")

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
        
        # 3. Define Headers for MAD Metrics
        mad_headers = OrderedDict()
        mad_headers['MAD of log ratios'] = {
            'title': 'MAD',
            'description': 'MAD QC: Median Absolute Deviation of log ratios',
            'format': '{:.4f}',
            'scale': 'RdYlGn',
            'hidden': False
        }
        mad_headers['Pearson correlation'] = {
            'title': 'Pearson',
            'description': 'MAD QC: Pearson correlation coefficient',
            'min': 0,
            'max': 1,
            'format': '{:.4f}',
            'scale': 'RdYlGn'
        }
        mad_headers['Spearman correlation'] = {
            'title': 'Spearman',
            'description': 'MAD QC: Spearman correlation coefficient',
            'min': 0,
            'max': 1,
            'format': '{:.4f}',
            'scale': 'RdYlGn'
        }

        # 4. Add to General Stats
        # We call this multiple times to merge data from different dictionaries
        self.general_stats_addcols(self.rnaseqc_data, rnaseqc_headers)
        self.general_stats_addcols(self.rsem_data, rsem_headers)
        self.general_stats_addcols(self.mad_data, mad_headers)

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
            plot=bargraph.plot(self.genetype_data, cats, pconfig)
        )