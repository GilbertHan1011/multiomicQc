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
            href='https://github.com/GilbertHan1011/multiomicQc',
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
        # Matches 'multiomics_report/gene_type_counts' from your execution_start
        for f in self.find_log_files('multiomics_report/gene_type_counts'):
            self.parse_genetype_json(f)

        # C. Parse RSEM Genes (JSON)
        # Matches 'multiomics/rsem' from your execution_start
        for f in self.find_log_files('multiomics_report/rsem'):
            self.parse_rsem_json(f)

        # D. Parse MAD files
        for f in self.find_log_files('multiomics_report/mad_rna'):
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
        """ Parses the nested Gene Type JSON file 
        
        Expected structure:
        {"gene_type_count": {"protein_coding": 1400, "lincRNA": 64, ...}}
        """
        try:
            data = json.loads(f['f'])
            # The file structure is {"gene_type_count": {...}}
            if "gene_type_count" in data:
                gene_counts = data["gene_type_count"]
                if isinstance(gene_counts, dict):
                    self.genetype_data[f['s_name']] = gene_counts
                else:
                    log.warning(f"gene_type_count is not a dict for sample {f['s_name']}, got type: {type(gene_counts)}")
            else:
                log.warning(f"Missing 'gene_type_count' key in JSON file for sample {f['s_name']}")
        except json.JSONDecodeError as e:
            log.warning(f"Failed to parse gene type JSON file {f.get('fn', 'unknown')}: {e}")
        except Exception as e:
            log.warning(f"Error parsing gene type file {f.get('fn', 'unknown')}: {e}")

    def _calculate_protein_coding_percentage(self):
        """
        Calculate the percentage of protein_coding genes for each sample.
        
        Returns:
            dict: {sample_name: {'protein_coding_percentage': value}}
        """
        protein_coding_data = {}
        
        for sample_name, gene_counts in self.genetype_data.items():
            if not isinstance(gene_counts, dict):
                continue
                
            # Get protein_coding count
            protein_coding_count = gene_counts.get('protein_coding', 0)
            
            # Calculate total count (sum of all gene types)
            total_count = sum(gene_counts.values())
            
            # Calculate percentage
            if total_count > 0:
                percentage = (protein_coding_count / total_count) * 100
                protein_coding_data[sample_name] = {
                    'protein_coding_percentage': percentage
                }
            else:
                log.warning(f"Sample {sample_name}: Total gene count is 0, cannot calculate protein_coding percentage")
        
        return protein_coding_data

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
            data = json.loads(f['f'])
            
            # Validation: specific key check
            if "MAD of log ratios" in data:
                # Use the clean sample name provided by MultiQC
                s_name = f['s_name']
                
                # Also check for status field to skip invalid files
                if "status" in data and data["status"] == "No pairs available for aggregation":
                    # Skip files with no pairs
                    log.warning(f"Skipping MAD QC file {f.get('fn', 'unknown')}: No pairs available")
                    return
                
                self.mad_data[s_name] = data
                self.add_data_source(f, section='mad_qc')
            else:
                log.warning(f"File {f.get('fn', 'unknown')} does not contain 'MAD of log ratios' key")
                
        except json.JSONDecodeError as e:
            log.warning(f"Failed to parse MAD QC JSON file {f.get('fn', 'unknown')}: {e}")
        except Exception as e:
            log.warning(f"Error parsing MAD QC file {f.get('fn', 'unknown')}: {e}")


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
        
        # 3. Calculate and add protein_coding percentage from gene type counts
        protein_coding_data = self._calculate_protein_coding_percentage()
        genetype_headers = OrderedDict()
        if protein_coding_data:
            genetype_headers['protein_coding_percentage'] = {
                'title': '% Protein Coding',
                'description': 'Percentage of protein coding genes among all detected gene types',
                'min': 0,
                'max': 100,
                'suffix': '%',
                'format': '{:,.1f}',
                'scale': 'YlGn'
            }
        
        # 4. Define Headers for MAD Metrics
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

        # 5. Add to General Stats
        # We call this multiple times to merge data from different dictionaries
        self.general_stats_addcols(self.rnaseqc_data, rnaseqc_headers)
        self.general_stats_addcols(self.rsem_data, rsem_headers)
        if protein_coding_data:
            self.general_stats_addcols(protein_coding_data, genetype_headers)
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
        cats = sorted(list(set(k for s in self.genetype_data.values() for k in s)))
        
        # Create the plot
        plot = bargraph.plot(self.genetype_data, cats, pconfig)
        
        if plot is None:
            log.warning("bargraph.plot() returned None - plot will not be displayed. Check data structure.")
            return
        
        self.add_section(
            name='Gene Types',
            anchor='my_rnaseq_genetypes',
            description='Counts of different gene biotypes (protein_coding, rRNA, etc.)',
            plot=plot
        )