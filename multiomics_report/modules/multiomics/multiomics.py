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
        self.peaks_stats_data = dict()  # Stores peaks stats data (peaks, frip, regulatory_fraction)
        self.prealign_data = dict()  # Stores prealign stats data (percent_filtered)
        self.atacseq_data = dict()  # Stores ATAC-seq/Samblaster data (n_tot, n_nondups, nrf)
        self.atacseq_tss_data = dict()  # Stores TSS coverage data (downsampled for plotting)
        self.replicate_correlations_data = dict()  # Stores replicate correlation data (pearson_p)
        self.reproducibility_qc_data = dict()  # Stores reproducibility QC data (rescue_ratio, self_consistency_ratio, N_optimal)
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
        
        # K. Parse Peaks Stats files (*stats.tsv)
        peaks_stats_files = list(self.find_log_files('multiomics_report/frip_atac'))
        log.debug(f"[DEBUG] Found {len(peaks_stats_files)} peaks stats files")
        for f in peaks_stats_files:
            self.parse_peaks_stats_tsv(f)
        
        # L. Parse Prealign Stats files (*prealign.stats.tsv)
        prealign_files = list(self.find_log_files('multiomics_report/prealign'))
        log.debug(f"[DEBUG] Found {len(prealign_files)} prealign stats files")
        for f in prealign_files:
            self.parse_prealign_stats_tsv(f)
        
        # M. Parse Samblaster files
        samblaster_files = list(self.find_log_files('multiomics_report/samblaster'))
        log.debug(f"[DEBUG] Found {len(samblaster_files)} samblaster files")
        for f in samblaster_files:
            self.parse_samblaster_log(f)
        
        # N. Parse TSS files
        tss_files = list(self.find_log_files('multiomics_report/tss'))
        log.debug(f"[DEBUG] Found {len(tss_files)} TSS files")
        for f in tss_files:
            my_sample_name = f['s_name'].replace('_TSS', '')
            tss_data, tss_max = self.parse_atacseq_tss(f['f'])
            self.atacseq_tss_data[f['s_name']] = tss_data
            # Store tss_max in atacseq_data with the cleaned sample name
            if my_sample_name not in self.atacseq_data:
                self.atacseq_data[my_sample_name] = {}
            self.atacseq_data[my_sample_name]['tss_max'] = tss_max
        
        # O. Parse Replicate Correlations files
        replicate_corr_files = list(self.find_log_files('multiomics_report/replicate_correlations'))
        log.debug(f"[DEBUG] Found {len(replicate_corr_files)} replicate correlation files")
        for f in replicate_corr_files:
            self.parse_replicate_correlations_tsv(f)
        
        # P. Parse Reproducibility QC files
        reproducibility_qc_files = list(self.find_log_files('multiomics_report/reproducibility_qc'))
        log.debug(f"[DEBUG] Found {len(reproducibility_qc_files)} reproducibility QC files")
        for f in reproducibility_qc_files:
            self.parse_reproducibility_qc_json(f)
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
        self.peaks_stats_data = self.ignore_samples(self.peaks_stats_data)
        self.prealign_data = self.ignore_samples(self.prealign_data)
        self.atacseq_data = self.ignore_samples(self.atacseq_data)
        self.atacseq_tss_data = self.ignore_samples(self.atacseq_tss_data)
        self.replicate_correlations_data = self.ignore_samples(self.replicate_correlations_data)
        self.reproducibility_qc_data = self.ignore_samples(self.reproducibility_qc_data)
        
        # If no data found at all, raise ModuleNoSamplesFound
        if (len(self.rnaseqc_data) == 0 and len(self.genetype_data) == 0 and 
            len(self.rsem_data) == 0 and len(self.mad_data) == 0 and
            len(self.coverage_data) == 0 and len(self.peak_count_data) == 0 and
            len(self.jaccard_data) == 0 and len(self.frip_data) == 0 and
            len(self.preseq_data) == 0 and len(self.bam_correlation_data) == 0 and
            len(self.peaks_stats_data) == 0 and len(self.prealign_data) == 0 and
            len(self.atacseq_data) == 0 and len(self.atacseq_tss_data) == 0 and
            len(self.replicate_correlations_data) == 0):
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

    def parse_peaks_stats_tsv(self, f):
        """ Parses the peaks stats TSV file 
        
        Expected format:
        peaks   1715
        frip    1
        regulatory_fraction     0.469257
        """
        try:
            parsed_data = {}
            for line in f['f'].splitlines():
                if not line.strip():
                    continue
                s = line.split('\t')
                if len(s) >= 2:
                    key = s[0].strip()
                    val = s[1].strip()
                    # Try to convert numbers to floats
                    try:
                        parsed_data[key] = float(val)
                    except ValueError:
                        parsed_data[key] = val  # Keep strings if not numeric
            
            # Only store if we have the expected keys
            if parsed_data and ('peaks' in parsed_data or 'frip' in parsed_data or 'regulatory_fraction' in parsed_data):
                # Rename keys to match expected format
                formatted_data = {}
                if 'peaks' in parsed_data:
                    formatted_data['peaks'] = int(parsed_data['peaks'])
                if 'frip' in parsed_data:
                    formatted_data['frip'] = parsed_data['frip']
                if 'regulatory_fraction' in parsed_data:
                    formatted_data['regulatory_fraction'] = parsed_data['regulatory_fraction']
                
                self.peaks_stats_data[f['s_name']] = formatted_data
                log.debug(f"Parsed peaks stats for '{f['s_name']}': {formatted_data}")
            else:
                log.warning(f"Peaks stats file {f.get('fn', 'unknown')} missing expected keys (peaks, frip, regulatory_fraction)")
                
        except Exception as e:
            log.warning(f"Error parsing peaks stats file {f.get('fn', 'unknown')}: {e}")

    def parse_prealign_stats_tsv(self, f):
        """ Parses the prealign stats TSV file 
        
        Expected format:
        prealignment    reads_before    reads_after     reads_filtered  percent_filtered
        chrM    2486    2160    326     0.131134
        total   2486    2160    326     0.131134
        """
        try:
            lines = f['f'].splitlines()
            if len(lines) < 2:
                log.warning(f"Prealign stats file {f.get('fn', 'unknown')} has insufficient lines")
                return
            
            # Parse header
            header = [col.strip() for col in lines[0].split('\t')]
            if 'percent_filtered' not in header:
                log.warning(f"Prealign stats file {f.get('fn', 'unknown')} missing 'percent_filtered' column")
                return
            
            percent_filtered_idx = header.index('percent_filtered')
            prealignment_idx = header.index('prealignment')
            
            # Find the row with "total" in the prealignment column
            percent_filtered = None
            for line_idx, line in enumerate(lines[1:], start=1):
                if not line.strip():
                    continue
                data_row = [val.strip() for val in line.split('\t')]
                if len(data_row) <= max(percent_filtered_idx, prealignment_idx):
                    continue
                
                # Check if this is the "total" row
                if data_row[prealignment_idx].lower() == 'total':
                    try:
                        percent_filtered = float(data_row[percent_filtered_idx])
                        # Convert decimal to percentage if value is < 1 (e.g., 0.131134 -> 13.11)
                        if percent_filtered < 1:
                            percent_filtered = percent_filtered * 100
                        break
                    except (ValueError, IndexError) as e:
                        log.warning(f"Prealign stats file {f.get('fn', 'unknown')}: Could not parse percent_filtered from total row: {e}")
                        continue
            
            if percent_filtered is not None:
                self.prealign_data[f['s_name']] = {'percent_filtered': percent_filtered}
                log.debug(f"Parsed prealign stats for '{f['s_name']}': percent_filtered={percent_filtered}")
            else:
                log.warning(f"Prealign stats file {f.get('fn', 'unknown')}: Could not find 'total' row or extract percent_filtered")
                
        except Exception as e:
            log.warning(f"Error parsing prealign stats file {f.get('fn', 'unknown')}: {e}")

    def parse_samblaster_log(self, f):
        """ Parses the Samblaster log file 
        
        Expected format:
        samblaster: Marked 3265 of 6480 (50.39%) read ids as duplicates using 2860k memory in 0.010S CPU seconds and 19S wall time.
        
        Extracts:
        - n_tot: total number of read ids
        - n_dups: number of duplicates
        - n_nondups: n_tot - n_dups
        - nrf: Non-Redundant Fraction = n_nondups / n_tot
        """
        import re
        try:
            sample_name = f['s_name']
            
            # Parse the log file
            for line in f['f'].splitlines():
                if 'Marked' in line and 'read ids as duplicates' in line:
                    # Example: "samblaster: Marked 3265 of 6480 (50.39%) read ids as duplicates..."
                    match = re.search(r'Marked\s+(\d+)\s+of\s+(\d+)', line)
                    if match:
                        n_dups = int(match.group(1))
                        n_tot = int(match.group(2))
                        n_nondups = n_tot - n_dups
                        
                        if sample_name not in self.atacseq_data:
                            self.atacseq_data[sample_name] = {}
                        
                        self.atacseq_data[sample_name]['n_tot'] = n_tot
                        self.atacseq_data[sample_name]['n_nondups'] = n_nondups
                        self.atacseq_data[sample_name]['n_dups'] = n_dups
                        
                        # Calculate NRF (Non-Redundant Fraction) = 1 - duplication_rate = n_nondups / n_tot
                        if n_tot > 0:
                            nrf = float(n_nondups) / float(n_tot)
                            self.atacseq_data[sample_name]['nrf'] = nrf
                        else:
                            self.atacseq_data[sample_name]['nrf'] = 0.0
                        
                        log.debug(f"Parsed samblaster data for '{sample_name}': n_tot={n_tot}, n_nondups={n_nondups}, nrf={nrf:.4f}")
                        break
            else:
                log.warning(f"Samblaster log file {f.get('fn', 'unknown')}: Could not find 'Marked ... read ids as duplicates' line")
                
        except Exception as e:
            log.warning(f"Error parsing samblaster log file {f.get('fn', 'unknown')}: {e}")

    def parse_atacseq_tss(self, file_content):
        """ Parses the TSS CSV file 
        
        Expected format:
        base,count
        -2001,0
        -2000,0
        ...
        
        Returns:
            tuple: (data, max_value)
                - data: OrderedDict mapping distance (int) to coverage (float), downsampled (every 10th point)
                - max_value: maximum coverage value
        """
        data = OrderedDict()
        count = 0
        max_value = 0.0
        
        try:
            for line in file_content.splitlines():
                if not line.strip():
                    continue
                    
                s = line.split(',')
                if len(s) < 2:
                    continue
                
                # Skip header line
                if s[0].strip() == 'base':
                    continue
                
                try:
                    distance = int(s[0].strip())
                    coverage = float(s[1].strip())
                    
                    # Track maximum coverage value
                    if coverage > max_value:
                        max_value = coverage
                    
                    count += 1
                    # Store every 10th data point for plotting (downsampling)
                    if count % 10 == 0:
                        data[distance] = coverage
                        
                except (ValueError, IndexError) as e:
                    log.warning(f"Could not parse TSS line: {line.strip()}, error: {e}")
                    continue
            
            log.debug(f"Parsed TSS data: {len(data)} points stored (downsampled), max_value={max_value}")
            return data, max_value
            
        except Exception as e:
            log.warning(f"Error parsing TSS file: {e}")
            return OrderedDict(), 0.0

    def parse_replicate_correlations_tsv(self, f):
        """ Parses the replicate correlations TSV file 
        
        Expected format:
        rep1    rep2    pearson_r       pearson_p       spearman_r      spearman_p
        test1-1 test1-2 0.9999999999999996      0.0     1.0     0.0
        
        Extracts pearson_p and stores it with the parent sample name (e.g., test1)
        """
        try:
            lines = f['f'].splitlines()
            if len(lines) < 2:
                log.warning(f"Replicate correlations file {f.get('fn', 'unknown')} has insufficient lines")
                return
            
            # Parse header
            header = [col.strip() for col in lines[0].split('\t')]
            if 'pearson_p' not in header:
                log.warning(f"Replicate correlations file {f.get('fn', 'unknown')} missing 'pearson_p' column")
                return
            
            pearson_p_idx = header.index('pearson_p')
            rep1_idx = header.index('rep1')
            
            # Parse data rows
            for line_idx, line in enumerate(lines[1:], start=1):
                if not line.strip():
                    continue
                    
                data_row = [val.strip() for val in line.split('\t')]
                if len(data_row) <= max(pearson_p_idx, rep1_idx):
                    continue
                
                try:
                    rep1 = data_row[rep1_idx]
                    pearson_p = float(data_row[pearson_p_idx])
                    
                    # Extract parent sample name (e.g., test1-1 -> test1)
                    # Split on hyphen and take the first part
                    parent_name = rep1.split('-')[0]
                    
                    # Store with parent sample name
                    if parent_name not in self.replicate_correlations_data:
                        self.replicate_correlations_data[parent_name] = {}
                    
                    self.replicate_correlations_data[parent_name]['pearson_p'] = pearson_p
                    log.debug(f"Parsed replicate correlations for parent '{parent_name}': pearson_p={pearson_p}")
                    
                except (ValueError, IndexError) as e:
                    log.warning(f"Replicate correlations file {f.get('fn', 'unknown')}: Could not parse row {line_idx}: {e}")
                    continue
                    
        except Exception as e:
            log.warning(f"Error parsing replicate correlations file {f.get('fn', 'unknown')}: {e}")

    def parse_reproducibility_qc_json(self, f):
        """ Parses the reproducibility QC JSON file 
        
        Expected format:
        {
          "prefix": "test1",
          "rescue_ratio": null,
          "self_consistency_ratio": 1.0,
          "N_optimal": 988
        }
        
        Extracts rescue_ratio, self_consistency_ratio, and N_optimal
        and stores them with the prefix (parent sample name)
        """
        try:
            data = json.loads(f['f'])
            
            # Extract prefix (parent sample name)
            if 'prefix' not in data:
                log.warning(f"Reproducibility QC file {f.get('fn', 'unknown')} missing 'prefix' key")
                return
            
            prefix = data['prefix']
            
            # Initialize storage for this sample
            if prefix not in self.reproducibility_qc_data:
                self.reproducibility_qc_data[prefix] = {}
            
            # Extract rescue_ratio (can be null)
            if 'rescue_ratio' in data:
                rescue_ratio = data['rescue_ratio']
                if rescue_ratio is not None:
                    try:
                        self.reproducibility_qc_data[prefix]['rescue_ratio'] = float(rescue_ratio)
                    except (ValueError, TypeError):
                        log.warning(f"Reproducibility QC file {f.get('fn', 'unknown')}: Could not convert rescue_ratio to float")
                # If null, we don't store it (or could store as None, but None values are typically skipped in MultiQC)
            
            # Extract self_consistency_ratio
            if 'self_consistency_ratio' in data:
                try:
                    self_consistency_ratio = float(data['self_consistency_ratio'])
                    self.reproducibility_qc_data[prefix]['self_consistency_ratio'] = self_consistency_ratio
                except (ValueError, TypeError) as e:
                    log.warning(f"Reproducibility QC file {f.get('fn', 'unknown')}: Could not parse self_consistency_ratio: {e}")
            
            # Extract N_optimal
            if 'N_optimal' in data:
                try:
                    n_optimal = int(data['N_optimal'])
                    self.reproducibility_qc_data[prefix]['N_optimal'] = n_optimal
                except (ValueError, TypeError) as e:
                    log.warning(f"Reproducibility QC file {f.get('fn', 'unknown')}: Could not parse N_optimal: {e}")
            
            log.debug(f"Parsed reproducibility QC for '{prefix}': {self.reproducibility_qc_data[prefix]}")
                
        except json.JSONDecodeError as e:
            log.warning(f"Failed to parse reproducibility QC JSON file {f.get('fn', 'unknown')}: {e}")
        except Exception as e:
            log.warning(f"Error parsing reproducibility QC file {f.get('fn', 'unknown')}: {e}")

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
        
        # 9. Define Headers for Peaks Stats
        peaks_stats_headers = OrderedDict()
        peaks_stats_headers['peaks'] = {
            'title': 'Peaks',
            'description': 'Number of peaks detected',
            'format': '{:,.0f}',
            'scale': 'Purples'
        }
        peaks_stats_headers['frip'] = {
            'title': 'FRIP',
            'description': 'Fraction of Reads in Peaks',
            'min': 0,
            'max': 1,
            'format': '{:.4f}',
            'scale': 'YlGn'
        }
        peaks_stats_headers['regulatory_fraction'] = {
            'title': 'Regulatory Fraction',
            'description': 'Fraction of peaks in regulatory regions',
            'min': 0,
            'max': 1,
            'format': '{:.2%}',
            'scale': 'Blues'
        }
        
        # 10. Define Headers for Prealign Stats
        prealign_headers = OrderedDict()
        prealign_headers['percent_filtered'] = {
            'title': 'Prealign Filtered %',
            'description': 'Percentage of reads filtered during prealignment',
            'min': 0,
            'max': 100,
            'suffix': '%',
            'format': '{:.2f}',
            'scale': 'Oranges'
        }
        
        # 11. Define Headers for Replicate Correlations
        replicate_corr_headers = OrderedDict()
        replicate_corr_headers['pearson_p'] = {
            'title': 'Peak Correlation',
            'description': 'Peak correlation p-value',
            'min': 0,
            'max': 1,
            'format': '{:.4f}',
            'scale': 'RdYlGn'
        }
        
        # 12. Define Headers for Reproducibility QC
        reproducibility_qc_headers = OrderedDict()
        reproducibility_qc_headers['rescue_ratio'] = {
            'title': 'Rescue Ratio',
            'description': 'IDR rescue ratio',
            'format': '{:.4f}',
            'scale': 'RdYlGn'
        }
        reproducibility_qc_headers['self_consistency_ratio'] = {
            'title': 'Self Consistency',
            'description': 'Self-consistency ratio',
            'min': 0,
            'max': 1,
            'format': '{:.4f}',
            'scale': 'RdYlGn'
        }
        reproducibility_qc_headers['N_optimal'] = {
            'title': 'N Optimal',
            'description': 'Number of optimal peaks',
            'format': '{:,.0f}',
            'scale': 'Blues'
        }
        
        # 13. Define Headers for Samblaster/ATAC-seq Stats
        samblaster_headers = OrderedDict()
        samblaster_headers['n_tot'] = {
            'title': 'Total\nAlignments',
            'description': 'Total number of alignments processed by Samblaster',
            'format': '{:,.0f}',
            'scale': 'Blues'
        }
        samblaster_headers['n_nondups'] = {
            'title': 'Non-Duplicate\nAlignments',
            'description': 'Number of non-duplicate alignments from Samblaster',
            'format': '{:,.0f}',
            'scale': 'Greens'
        }
        samblaster_headers['nrf'] = {
            'title': 'NRF',
            'description': 'Non-Redundant Fraction (1 - duplication rate)',
            'min': 0.0,
            'max': 1.0,
            'format': '{:.4f}',
            'scale': 'Greens'
        }
        samblaster_headers['tss_max'] = {
            'title': 'TSS Max',
            'description': 'Maximum TSS enrichment value',
            'format': '{:.2f}',
            'scale': 'RdYlGn',
            'min': 0
        }

        # 14. Add to General Stats
        # We call this multiple times to merge data from different dictionaries
        self.general_stats_addcols(self.rnaseqc_data, rnaseqc_headers)
        self.general_stats_addcols(self.rsem_data, rsem_headers)
        self.general_stats_addcols(self.coverage_data, coverage_headers)
        self.general_stats_addcols(self.peak_count_data, peak_count_headers)
        self.general_stats_addcols(self.jaccard_data, jaccard_headers)
        self.general_stats_addcols(self.frip_data, frip_headers)
        self.general_stats_addcols(self.preseq_data, preseq_headers)
        self.general_stats_addcols(self.bam_correlation_data, bam_correlation_headers)
        self.general_stats_addcols(self.peaks_stats_data, peaks_stats_headers)
        self.general_stats_addcols(self.prealign_data, prealign_headers)
        self.general_stats_addcols(self.replicate_correlations_data, replicate_corr_headers)
        self.general_stats_addcols(self.reproducibility_qc_data, reproducibility_qc_headers)
        
        # Prepare samblaster data with proper type conversion
        samblaster_data = {}
        for sample_name in self.atacseq_data:
            samblaster_data[sample_name] = {}
            
            if 'n_tot' in self.atacseq_data[sample_name]:
                try:
                    value = int(self.atacseq_data[sample_name]['n_tot'])
                except (ValueError, TypeError):
                    value = None
                samblaster_data[sample_name]['n_tot'] = value
            
            if 'n_nondups' in self.atacseq_data[sample_name]:
                try:
                    value = int(self.atacseq_data[sample_name]['n_nondups'])
                except (ValueError, TypeError):
                    value = None
                samblaster_data[sample_name]['n_nondups'] = value
            
            if 'nrf' in self.atacseq_data[sample_name]:
                try:
                    value = float(self.atacseq_data[sample_name]['nrf'])
                except (ValueError, TypeError):
                    value = None
                samblaster_data[sample_name]['nrf'] = value
            
            if 'tss_max' in self.atacseq_data[sample_name]:
                try:
                    value = float(self.atacseq_data[sample_name]['tss_max'])
                except (ValueError, TypeError):
                    value = None
                samblaster_data[sample_name]['tss_max'] = value
        
        self.general_stats_addcols(samblaster_data, samblaster_headers)
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