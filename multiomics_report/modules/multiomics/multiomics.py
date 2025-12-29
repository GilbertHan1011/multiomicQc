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
        self.coverage_data = dict()   # Stores coverage.tsv data
        self.peak_count_data = dict() # Stores peak count data
        self.jaccard_data = dict()    # Stores jaccard data
        self.frip_data = dict()       # Stores FRIP data
        self.preseq_data = dict()     # Stores preseq data
        self.bam_correlation_data = dict()  # Stores bam_correlation data
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
        
        # E. Parse Coverage files
        coverage_files = list(self.find_log_files('multiomics_report/peak_coverage'))
        log.debug(f"[DEBUG] Found {len(coverage_files)} coverage files")
        for f in coverage_files:
            self.parse_coverage_tsv(f)
        
        # F. Parse Peak Count files
        peak_count_files = list(self.find_log_files('multiomics_report/peak_count'))
        log.debug(f"[DEBUG] Found {len(peak_count_files)} peak count files")
        for f in peak_count_files:
            self.parse_peak_count_txt(f)
        
        # G. Parse Jaccard files
        jaccard_files = list(self.find_log_files('multiomics_report/jaccard'))
        log.debug(f"[DEBUG] Found {len(jaccard_files)} jaccard files")
        for f in jaccard_files:
            self.parse_jaccard_txt(f)
        
        # H. Parse FRIP files
        frip_files = list(self.find_log_files('multiomics_report/frip'))
        log.debug(f"[DEBUG] Found {len(frip_files)} FRIP files")
        for f in frip_files:
            self.parse_frip_tsv(f)
        
        # I. Parse Preseq files
        for f in self.find_log_files('multiomics_report/preseq'):
            self.parse_preseq_txt(f)
        
        # J. Parse BAM Correlation files
        bam_correlation_files = list(self.find_log_files('multiomics_report/bam_correlation'))
        log.debug(f"[DEBUG] Found {len(bam_correlation_files)} bam_correlation files")
        for f in bam_correlation_files:
            self.parse_bam_correlation_tsv(f)
        # -----------------------------------------------------------
        # 4. FILTERING & EXIT
        # -----------------------------------------------------------
        # Filter out samples ignored by the user
        self.rnaseqc_data = self.ignore_samples(self.rnaseqc_data)
        self.genetype_data = self.ignore_samples(self.genetype_data)
        self.rsem_data = self.ignore_samples(self.rsem_data)
        self.mad_data = self.ignore_samples(self.mad_data)
        self.coverage_data = self.ignore_samples(self.coverage_data)
        self.peak_count_data = self.ignore_samples(self.peak_count_data)
        self.jaccard_data = self.ignore_samples(self.jaccard_data)
        self.frip_data = self.ignore_samples(self.frip_data)
        self.preseq_data = self.ignore_samples(self.preseq_data)
        self.bam_correlation_data = self.ignore_samples(self.bam_correlation_data)
        
        # If no data found at all, raise ModuleNoSamplesFound
        if (len(self.rnaseqc_data) == 0 and len(self.genetype_data) == 0 and 
            len(self.rsem_data) == 0 and len(self.mad_data) == 0 and
            len(self.coverage_data) == 0 and len(self.peak_count_data) == 0 and
            len(self.jaccard_data) == 0 and len(self.frip_data) == 0 and
            len(self.preseq_data) == 0 and len(self.bam_correlation_data) == 0):
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

    def parse_coverage_tsv(self, f):
        """ Parses the coverage TSV file 
        
        Expected format:
        Sample	Coverage_Percent	Covered_Bases	Total_Genome_Size
        P6-1_test	.035600	1100659	3088286401
        """
        try:
            lines = f['f'].splitlines()
            if len(lines) < 2:
                log.warning(f"Coverage file {f.get('fn', 'unknown')} has insufficient lines")
                return
            
            # Parse header
            header = lines[0].split('\t')
            if 'Coverage_Percent' not in header:
                log.warning(f"Coverage file {f.get('fn', 'unknown')} missing 'Coverage_Percent' column")
                return
            
            # Parse data row (assuming one sample per file)
            data_row = lines[1].split('\t')
            if len(data_row) != len(header):
                log.warning(f"Coverage file {f.get('fn', 'unknown')} has mismatched columns")
                return
            
            # Create dict from header and data
            parsed_data = {}
            for i, col in enumerate(header):
                col = col.strip()
                try:
                    # Try to convert to float, keep as string if it fails
                    parsed_data[col] = float(data_row[i].strip())
                except ValueError:
                    parsed_data[col] = data_row[i].strip()
            
            # Store with sample name
            self.coverage_data[f['s_name']] = parsed_data
            
        except Exception as e:
            log.warning(f"Error parsing coverage file {f.get('fn', 'unknown')}: {e}")

    def parse_peak_count_txt(self, f):
        """ Parses the peak count text file 
        
        Expected format:
        P6-2_test 3984
        """
        try:
            content = f['f'].strip()
            # Try to extract number (could be on same line as sample name or separate)
            parts = content.split()
            
            # Find the numeric value
            peak_count = None
            for part in parts:
                try:
                    peak_count = int(float(part))
                    break
                except ValueError:
                    continue
            
            if peak_count is not None:
                self.peak_count_data[f['s_name']] = {'peak_count': peak_count}
            else:
                log.warning(f"Could not extract peak count from file {f.get('fn', 'unknown')}")
                
        except Exception as e:
            log.warning(f"Error parsing peak count file {f.get('fn', 'unknown')}: {e}")

    def parse_jaccard_txt(self, f):
        """ Parses the jaccard text file 
        
        Expected format:
        Sample1	Sample2	Jaccard_Similarity	Intersection	Union_Intersection
        P6-1_test	P6-2_test	1	1100659	1100659
        (may have multiple rows - will calculate mean of numeric columns)
        """
        try:
            lines = f['f'].splitlines()
            if len(lines) < 2:
                log.warning(f"Jaccard file {f.get('fn', 'unknown')} has insufficient lines")
                return
            
            # Parse header
            header = [col.strip() for col in lines[0].split('\t')]
            if 'Jaccard_Similarity' not in header:
                log.warning(f"Jaccard file {f.get('fn', 'unknown')} missing 'Jaccard_Similarity' column")
                return
            
            # Identify numeric columns (exclude Sample1, Sample2)
            numeric_columns = ['Jaccard_Similarity', 'Intersection', 'Union_Intersection']
            numeric_col_indices = [i for i, col in enumerate(header) if col in numeric_columns]
            
            # Collect all data rows
            all_rows = []
            for line_idx, line in enumerate(lines[1:], start=1):
                if not line.strip():
                    continue
                data_row = [val.strip() for val in line.split('\t')]
                if len(data_row) != len(header):
                    log.warning(f"Jaccard file {f.get('fn', 'unknown')} row {line_idx} has mismatched columns (expected {len(header)}, got {len(data_row)})")
                    continue
                all_rows.append(data_row)
            
            if not all_rows:
                log.warning(f"Jaccard file {f.get('fn', 'unknown')} has no valid data rows")
                return
            
            # Calculate mean for numeric columns
            parsed_data = {}
            
            # For numeric columns, calculate mean
            for col_idx in numeric_col_indices:
                col_name = header[col_idx]
                numeric_values = []
                for row in all_rows:
                    try:
                        numeric_values.append(float(row[col_idx]))
                    except (ValueError, IndexError):
                        log.warning(f"Jaccard file {f.get('fn', 'unknown')}: Could not parse numeric value for column '{col_name}' in row")
                        continue
                
                if numeric_values:
                    parsed_data[col_name] = sum(numeric_values) / len(numeric_values)
                else:
                    log.warning(f"Jaccard file {f.get('fn', 'unknown')}: No valid numeric values for column '{col_name}'")
            
            # For non-numeric columns (Sample1, Sample2), use first row value
            for col_idx, col_name in enumerate(header):
                if col_name not in numeric_columns:
                    if all_rows:
                        parsed_data[col_name] = all_rows[0][col_idx]
            
            # Store with sample name
            s_name = f['s_name']
            self.jaccard_data[s_name] = parsed_data
            if len(all_rows) > 1:
                log.debug(f"Plugin: Parsed jaccard data for '{s_name}' from {len(all_rows)} rows (mean calculated): {parsed_data}")
            else:
                log.debug(f"Plugin: Parsed jaccard data for '{s_name}': {parsed_data}")
            
        except Exception as e:
            log.warning(f"Error parsing jaccard file {f.get('fn', 'unknown')}: {e}")

    def parse_frip_tsv(self, f):
        """ Parses the FRIP TSV file 
        
        Expected format:
        file	featureType	percent	featureReadCount	totalReadCount
        /path/to/file.bam	frip	4.72	45984	975200
        """
        try:
            lines = f['f'].splitlines()
            if len(lines) < 2:
                log.warning(f"FRIP file {f.get('fn', 'unknown')} has insufficient lines")
                return
            
            # Parse header
            header = lines[0].split('\t')
            if 'percent' not in header:
                log.warning(f"FRIP file {f.get('fn', 'unknown')} missing 'percent' column")
                return
            
            # Parse data rows (could be multiple samples in one file)
            for line_idx, line in enumerate(lines[1:], start=1):
                if not line.strip():
                    continue
                    
                data_row = line.split('\t')
                if len(data_row) != len(header):
                    log.warning(f"FRIP file {f.get('fn', 'unknown')} row {line_idx} has mismatched columns")
                    continue
                
                # Create dict from header and data
                parsed_data = {}
                for i, col in enumerate(header):
                    col = col.strip()
                    try:
                        parsed_data[col] = float(data_row[i].strip())
                    except ValueError:
                        parsed_data[col] = data_row[i].strip()
                
                # Check if featureType is 'frip'
                if parsed_data.get('featureType', '').lower() == 'frip':
                    # Use sample name from file, or extract from file path if needed
                    s_name = f['s_name']
                    self.frip_data[s_name] = {'percent': parsed_data['percent']}
            
        except Exception as e:
            log.warning(f"Error parsing FRIP file {f.get('fn', 'unknown')}: {e}")

    def parse_preseq_txt(self, f):
        """ Parses the preseq lcextrap output file 
        
        Extracts the 200th row value (expected distinct reads at ~200M depth).
        This represents the library complexity at high sequencing depth.
        
        Expected format:
        TOTAL_READS	EXPECTED_DISTINCT
        1000	950
        ...
        (200th data row)
        """
        try:
            lines = f['f'].splitlines()
            if len(lines) < 201:  # Need header + 200 data rows
                log.warning(f"Preseq file {f.get('fn', 'unknown')} has insufficient lines (need at least 201, got {len(lines)})")
                return
            
            # Get the 200th row (line 201 including header, index 200)
            row_200 = lines[200].strip().split('\t')
            if len(row_200) < 2:
                log.warning(f"Preseq file {f.get('fn', 'unknown')} row 200 has insufficient columns")
                return
            
            try:
                # Column 0: TOTAL_READS, Column 1: EXPECTED_DISTINCT
                expected_distinct = float(row_200[1])
                s_name = f['s_name']
                self.preseq_data[s_name] = {'preseq_200M': expected_distinct}
                log.debug(f"Added preseq 200th row for sample '{s_name}': {expected_distinct:,.0f}")
            except (ValueError, IndexError) as e:
                log.warning(f"Could not parse preseq value from {f.get('fn', 'unknown')}: {e}")
                
        except Exception as e:
            log.warning(f"Error parsing preseq file {f.get('fn', 'unknown')}: {e}")

    def parse_bam_correlation_tsv(self, f):
        """ Parses the BAM correlation TSV file 
        
        Expected format:
        Sample	BAM_Correlation	Correlation_Quality	N_Replicates
        P6	1.0000	excellent	2
        """
        try:
            lines = f['f'].splitlines()
            if len(lines) < 2:
                log.warning(f"BAM correlation file {f.get('fn', 'unknown')} has insufficient lines")
                return
            
            # Parse header
            header = [col.strip() for col in lines[0].split('\t')]
            if 'BAM_Correlation' not in header:
                log.warning(f"BAM correlation file {f.get('fn', 'unknown')} missing 'BAM_Correlation' column")
                return
            
            # Parse data rows (could be multiple samples in one file)
            for line_idx, line in enumerate(lines[1:], start=1):
                if not line.strip():
                    continue
                    
                data_row = [val.strip() for val in line.split('\t')]
                if len(data_row) != len(header):
                    log.warning(f"BAM correlation file {f.get('fn', 'unknown')} row {line_idx} has mismatched columns (expected {len(header)}, got {len(data_row)})")
                    continue
                
                # Create dict from header and data
                parsed_data = {}
                for i, col in enumerate(header):
                    col = col.strip()
                    try:
                        # Try to convert to float for numeric columns
                        parsed_data[col] = float(data_row[i].strip())
                    except ValueError:
                        # Keep as string for non-numeric columns (e.g., Correlation_Quality)
                        parsed_data[col] = data_row[i].strip()
                
                # Extract sample name from the data (use Sample column if available, otherwise use file's s_name)
                s_name = parsed_data.get('Sample', f['s_name'])
                
                # Store only BAM_Correlation for general stats
                self.bam_correlation_data[s_name] = {'BAM_Correlation': parsed_data.get('BAM_Correlation')}
                log.debug(f"Parsed BAM correlation data for '{s_name}': BAM_Correlation={parsed_data.get('BAM_Correlation')}")
            
        except Exception as e:
            log.warning(f"Error parsing BAM correlation file {f.get('fn', 'unknown')}: {e}")

    # ===============================================================
    # REPORT WRITING
    # ===============================================================

    def write_general_stats(self):
        """ Adds columns to the main top table in MultiQC """
        
        # 1. Define Headers for RNA-SeQC Metrics
        rnaseqc_headers = OrderedDict()
        rnaseqc_headers['Exonic Rate'] = {
            'title': 'Exonic',
            'description': 'RNA-SeQC: Exonic Rate',
            'max': 1, 'min': 0,
            'scale': 'Greens',
            'format': '{:,.1%}'
        }
        rnaseqc_headers['Intergenic Rate'] = {
            'title': 'Intergenic',
            'description': 'RNA-SeQC: Intergenic Rate',
            'max': 1, 'min': 0,
            'scale': 'Oranges',
            'format': '{:,.1%}'
        }
        rnaseqc_headers['rRNA Rate'] = {
            'title': 'rRNA',
            'description': 'RNA-SeQC: rRNA Rate',
            'max': 1, 'min': 0,
            'scale': 'Reds',
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
        
        # 3. Define Headers for Coverage Metrics
        coverage_headers = OrderedDict()
        coverage_headers['Coverage_Percent'] = {
            'title': 'Coverage %',
            'description': 'Peak Coverage: Percentage of genome covered',
            'min': 0,
            'max': 100,
            'suffix': '%',
            'format': '{:.4f}',
            'scale': 'Blues'
        }
        
        # 4. Define Headers for Peak Count
        peak_count_headers = OrderedDict()
        peak_count_headers['peak_count'] = {
            'title': 'Peak Count',
            'description': 'Number of peaks detected',
            'format': '{:,.0f}',
            'scale': 'Purples'
        }
        
        # 5. Define Headers for Jaccard
        jaccard_headers = OrderedDict()
        jaccard_headers['Jaccard_Similarity'] = {
            'title': 'Jaccard',
            'description': 'Jaccard similarity coefficient',
            'min': 0,
            'max': 1,
            'format': '{:.4f}',
            'scale': 'RdYlGn'
        }
        # 6. Define Headers for FRIP
        frip_headers = OrderedDict()
        frip_headers['percent'] = {
            'title': 'FRIP %',
            'description': 'Fraction of Reads in Peaks (FRIP) percentage',
            'min': 0,
            'max': 100,
            'suffix': '%',
            'format': '{:.2f}',
            'scale': 'YlGn'
        }
        
        # 7. Define Headers for Preseq
        preseq_headers = OrderedDict()
        preseq_headers['preseq_200M'] = {
            'title': 'Preseq 200M',
            'description': 'Expected distinct reads at 200M depth (library complexity)',
            'format': '{:,.0f}',
            'scale': 'RdYlGn',
            'min': 0
        }
        
        # 8. Define Headers for BAM Correlation
        bam_correlation_headers = OrderedDict()
        bam_correlation_headers['BAM_Correlation'] = {
            'title': 'BAM Correlation',
            'description': 'BAM file correlation coefficient between replicates',
            'min': 0,
            'max': 1,
            'format': '{:.4f}',
            'scale': 'RdYlGn'
        }

        # 9. Add to General Stats
        # We call this multiple times to merge data from different dictionaries
        self.general_stats_addcols(self.rnaseqc_data, rnaseqc_headers)
        self.general_stats_addcols(self.rsem_data, rsem_headers)
        self.general_stats_addcols(self.coverage_data, coverage_headers)
        self.general_stats_addcols(self.peak_count_data, peak_count_headers)
        self.general_stats_addcols(self.jaccard_data, jaccard_headers)
        self.general_stats_addcols(self.frip_data, frip_headers)
        self.general_stats_addcols(self.preseq_data, preseq_headers)
        self.general_stats_addcols(self.bam_correlation_data, bam_correlation_headers)
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