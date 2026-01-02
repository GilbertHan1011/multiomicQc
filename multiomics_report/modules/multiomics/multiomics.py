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
        self.hic_mapstat_data = dict()  # Stores Hi-C BAM mapping statistics (mapped_rate, etc.)
        self.hic_pairstat_data = dict()  # Stores Hi-C pair alignment statistics
        self.hic_rsstat_data = dict()  # Stores Hi-C processing statistics (RSstat)
        self.hic_dedup_stats_data = dict()  # Stores Hi-C deduplication statistics
        self.hic_dist_contact_data = dict()  # Stores Hi-C distance vs contact statistics (slope, mse)
        self.hic_loop_counts_data = dict()  # Stores Hi-C loop counts statistics (total loop count)
        self.hic_library_complexity_data = dict()  # Stores Hi-C library complexity statistics (C)
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
        
        # Q. Parse Hi-C BAM mapping statistics (mapstat files)
        hic_mapstat_files = list(self.find_log_files('multiomics_report/hic_mapstat'))
        log.debug(f"[DEBUG] Found {len(hic_mapstat_files)} Hi-C mapstat files")
        for f in hic_mapstat_files:
            self.parse_hic_mapstat(f)
        
        # R. Parse Hi-C pair alignment statistics (pairstat files)
        hic_pairstat_files = list(self.find_log_files('multiomics_report/hic_pairstat'))
        log.debug(f"[DEBUG] Found {len(hic_pairstat_files)} Hi-C pairstat files")
        for f in hic_pairstat_files:
            self.parse_hic_pairstat(f)
        
        # S. Parse Hi-C processing statistics (RSstat files)
        hic_rsstat_files = list(self.find_log_files('multiomics_report/RSstat'))
        log.debug(f"[DEBUG] Found {len(hic_rsstat_files)} Hi-C RSstat files")
        for f in hic_rsstat_files:
            self.parse_hic_rsstat(f)
        
        # T. Parse Hi-C deduplication statistics (dedup.stats files)
        hic_dedup_files = list(self.find_log_files('multiomics_report/dedup_stats'))
        log.debug(f"[DEBUG] Found {len(hic_dedup_files)} Hi-C dedup stats files")
        for f in hic_dedup_files:
            self.parse_hic_dedup_stats(f)
        
        # U. Parse Hi-C distance vs contact statistics (loglog_fits.csv files)
        hic_dist_contact_files = list(self.find_log_files('multiomics_report/hic_dist_contact'))
        log.debug(f"[DEBUG] Found {len(hic_dist_contact_files)} Hi-C distance vs contact files")
        for f in hic_dist_contact_files:
            self.parse_hic_dist_contact_csv(f)
        
        # V. Parse Hi-C loop counts statistics (loop_counts.tsv files)
        hic_loop_counts_files = list(self.find_log_files('multiomics_report/hic_loop_counts'))
        log.debug(f"[DEBUG] Found {len(hic_loop_counts_files)} Hi-C loop counts files")
        for f in hic_loop_counts_files:
            self.parse_hic_loop_counts_tsv(f)
        
        # W. Parse Hi-C library complexity statistics (complexity*.tsv files)
        hic_complexity_files = list(self.find_log_files('multiomics_report/hic_library_complexity'))
        log.debug(f"[DEBUG] Found {len(hic_complexity_files)} Hi-C library complexity files")
        for f in hic_complexity_files:
            self.parse_hic_library_complexity_tsv(f)
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
        self.hic_mapstat_data = self.ignore_samples(self.hic_mapstat_data)
        self.hic_pairstat_data = self.ignore_samples(self.hic_pairstat_data)
        self.hic_rsstat_data = self.ignore_samples(self.hic_rsstat_data)
        self.hic_dedup_stats_data = self.ignore_samples(self.hic_dedup_stats_data)
        self.hic_dist_contact_data = self.ignore_samples(self.hic_dist_contact_data)
        self.hic_loop_counts_data = self.ignore_samples(self.hic_loop_counts_data)
        self.hic_library_complexity_data = self.ignore_samples(self.hic_library_complexity_data)
        
        # If no data found at all, raise ModuleNoSamplesFound
        if (len(self.rnaseqc_data) == 0 and len(self.genetype_data) == 0 and 
            len(self.rsem_data) == 0 and len(self.mad_data) == 0 and
            len(self.coverage_data) == 0 and len(self.peak_count_data) == 0 and
            len(self.jaccard_data) == 0 and len(self.frip_data) == 0 and
            len(self.preseq_data) == 0 and len(self.bam_correlation_data) == 0 and
            len(self.peaks_stats_data) == 0 and len(self.prealign_data) == 0 and
            len(self.atacseq_data) == 0 and len(self.atacseq_tss_data) == 0 and
            len(self.replicate_correlations_data) == 0 and len(self.hic_mapstat_data) == 0 and
            len(self.hic_pairstat_data) == 0 and len(self.hic_rsstat_data) == 0 and
            len(self.hic_dedup_stats_data) == 0 and len(self.hic_dist_contact_data) == 0 and
            len(self.hic_loop_counts_data) == 0 and len(self.hic_library_complexity_data) == 0):
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

    def parse_hic_mapstat(self, f):
        """ Parses the Hi-C BAM mapping statistics file (mapstat)
        
        Expected format:
        ## /path/to/file
        total_1 24859
        mapped_1        24155
        global_1        9174
        local_1 14981
        
        Calculates: mapped_rate = mapped_1 / total_1
        """
        try:
            parsed_data = {}
            total_1 = None
            mapped_1 = None
            
            for line in f['f'].splitlines():
                if not line.strip() or line.strip().startswith('##'):
                    continue
                
                parts = line.split()
                if len(parts) >= 2:
                    key = parts[0].strip()
                    try:
                        value = int(parts[1].strip())
                        parsed_data[key] = value
                        
                        if key == 'total_1':
                            total_1 = value
                        elif key == 'mapped_1':
                            mapped_1 = value
                    except ValueError:
                        # Skip non-numeric values
                        continue
            
            # Calculate mapped_rate
            if total_1 is not None and total_1 > 0 and mapped_1 is not None:
                mapped_rate = float(mapped_1) / float(total_1)
                parsed_data['mapped_rate'] = mapped_rate
                log.debug(f"Calculated mapped_rate for '{f['s_name']}': {mapped_rate:.4f} ({mapped_1}/{total_1})")
            
            if parsed_data:
                self.hic_mapstat_data[f['s_name']] = parsed_data
                log.debug(f"Parsed Hi-C mapstat for '{f['s_name']}': {parsed_data}")
            else:
                log.warning(f"Hi-C mapstat file {f.get('fn', 'unknown')} has no valid data")
                
        except Exception as e:
            log.warning(f"Error parsing Hi-C mapstat file {f.get('fn', 'unknown')}: {e}")

    def parse_hic_pairstat(self, f):
        """ Parses the Hi-C pair alignment statistics file (pairstat)
        
        Expected format:
        Total_pairs_processed   24859   100.0
        Unmapped_pairs  197     0.792
        Low_qual_pairs  9787    39.37
        Unique_paired_alignments        13801   55.517
        Multiple_pairs_alignments       0       0.0
        Pairs_with_singleton    1074    4.32
        Low_qual_singleton      0       0.0
        Unique_singleton_alignments     0       0.0
        Multiple_singleton_alignments   0       0.0
        Reported_pairs  13801   55.517
        """
        try:
            parsed_data = {}
            
            for line in f['f'].splitlines():
                if not line.strip():
                    continue
                
                parts = line.split('\t')
                if len(parts) >= 2:
                    key = parts[0].strip()
                    try:
                        # Try to get the count (second column) and percentage (third column if available)
                        count = int(float(parts[1].strip()))
                        parsed_data[key] = count
                        
                        # Also store percentage if available
                        if len(parts) >= 3:
                            try:
                                percentage = float(parts[2].strip())
                                parsed_data[f"{key}_percent"] = percentage
                            except ValueError:
                                pass
                    except (ValueError, IndexError):
                        # Skip if not numeric
                        continue
            
            # Calculate Mapped_ratio = 1 - (Unmapped_pairs / Total_pairs_processed)
            if 'Unmapped_pairs' in parsed_data and 'Total_pairs_processed' in parsed_data:
                total_pairs = parsed_data['Total_pairs_processed']
                if total_pairs > 0:
                    unmapped = parsed_data['Unmapped_pairs']
                    mapped_ratio = 1.0 - (float(unmapped) / float(total_pairs))
                    parsed_data['Mapped_ratio'] = mapped_ratio
                    log.debug(f"Calculated Mapped_ratio for '{f['s_name']}': {mapped_ratio:.4f}")
            
            # Calculate unique_ratio = Unique_paired_alignments / Total_pairs_processed
            if 'Unique_paired_alignments' in parsed_data and 'Total_pairs_processed' in parsed_data:
                total_pairs = parsed_data['Total_pairs_processed']
                if total_pairs > 0:
                    unique_pairs = parsed_data['Unique_paired_alignments']
                    unique_ratio = float(unique_pairs) / float(total_pairs)
                    parsed_data['unique_ratio'] = unique_ratio
                    log.debug(f"Calculated unique_ratio for '{f['s_name']}': {unique_ratio:.4f}")
            
            if parsed_data:
                self.hic_pairstat_data[f['s_name']] = parsed_data
                log.debug(f"Parsed Hi-C pairstat for '{f['s_name']}': {parsed_data}")
            else:
                log.warning(f"Hi-C pairstat file {f.get('fn', 'unknown')} has no valid data")
                
        except Exception as e:
            log.warning(f"Error parsing Hi-C pairstat file {f.get('fn', 'unknown')}: {e}")

    def parse_hic_rsstat(self, f):
        """ Parses the Hi-C processing statistics file (RSstat)
        
        Expected format:
        ## Hi-C processing
        Valid_interaction_pairs 13247
        Valid_interaction_pairs_FF      3374
        Valid_interaction_pairs_RR      3293
        Valid_interaction_pairs_RF      3228
        Valid_interaction_pairs_FR      3352
        Dangling_end_pairs      443
        Religation_pairs        103
        Self_Cycle_pairs        8
        Single-end_pairs        0
        Filtered_pairs  0
        Dumped_pairs    0
        """
        try:
            parsed_data = {}
            
            for line in f['f'].splitlines():
                if not line.strip() or line.strip().startswith('##'):
                    continue
                
                parts = line.split()
                if len(parts) >= 2:
                    key = parts[0].strip()
                    try:
                        value = int(parts[1].strip())
                        parsed_data[key] = value
                    except ValueError:
                        # Skip non-numeric values
                        continue
            
            if parsed_data:
                self.hic_rsstat_data[f['s_name']] = parsed_data
                log.debug(f"Parsed Hi-C RSstat for '{f['s_name']}': {parsed_data}")
            else:
                log.warning(f"Hi-C RSstat file {f.get('fn', 'unknown')} has no valid data")
                
        except Exception as e:
            log.warning(f"Error parsing Hi-C RSstat file {f.get('fn', 'unknown')}: {e}")

    def parse_hic_dedup_stats(self, f):
        """ Parses the Hi-C deduplication statistics file (dedup.stats)
        
        Expected format:
        total   13801
        total_unmapped  0
        total_single_sided_mapped       0
        total_mapped    13801
        total_dups      1475
        total_nodups    12326
        cis     9056
        trans   3270
        cis_2kb+        5000
        pair_types/UU   12326
        pair_types/DD   1475
        summary/frac_dups   0.1069
        """
        try:
            parsed_data = {}
            
            for line in f['f'].splitlines():
                if not line.strip():
                    continue
                
                parts = line.split()
                if len(parts) >= 2:
                    key = parts[0].strip()
                    try:
                        value = int(float(parts[1].strip()))
                        parsed_data[key] = value
                    except ValueError:
                        # Try as float for frac_dups
                        try:
                            value = float(parts[1].strip())
                            parsed_data[key] = value
                        except ValueError:
                            # Skip non-numeric values
                            continue
            
            # Handle cis_2kb+ or cis_2kb field name variations
            cis_2kb_value = None
            for key in ['cis_2kb+', 'cis_2kb', 'cis_2kb_plus']:
                if key in parsed_data:
                    cis_2kb_value = parsed_data[key]
                    # Normalize to cis_2kb for easier access
                    if key != 'cis_2kb':
                        parsed_data['cis_2kb'] = cis_2kb_value
                    break
            
            # Normalize summary/frac_dups to frac_dups
            if 'summary/frac_dups' in parsed_data:
                parsed_data['frac_dups'] = parsed_data['summary/frac_dups']
            
            # Calculate frac_trans = trans / total
            if 'trans' in parsed_data and 'total' in parsed_data:
                total = parsed_data['total']
                if total > 0:
                    frac_trans = float(parsed_data['trans']) / float(total)
                    parsed_data['frac_trans'] = frac_trans
                    log.debug(f"Calculated frac_trans for '{f['s_name']}': {frac_trans:.4f}")
            
            # Calculate cis_trans_ratio = cis / trans
            if 'cis' in parsed_data and 'trans' in parsed_data:
                trans = parsed_data['trans']
                if trans > 0:
                    cis_trans_ratio = float(parsed_data['cis']) / float(trans)
                    parsed_data['cis_trans_ratio'] = cis_trans_ratio
                    log.debug(f"Calculated cis_trans_ratio for '{f['s_name']}': {cis_trans_ratio:.4f}")
            
            if parsed_data:
                self.hic_dedup_stats_data[f['s_name']] = parsed_data
                log.debug(f"Parsed Hi-C dedup stats for '{f['s_name']}': {parsed_data}")
            else:
                log.warning(f"Hi-C dedup stats file {f.get('fn', 'unknown')} has no valid data")
                
        except Exception as e:
            log.warning(f"Error parsing Hi-C dedup stats file {f.get('fn', 'unknown')}: {e}")

    def parse_hic_dist_contact_csv(self, f):
        """ Parses the Hi-C distance vs contact statistics CSV file (loglog_fits.csv)
        
        Expected format:
        region,slope,mse,n_points
        ALL_REGIONS,,,0
        
        Extracts slope and mse from the ALL_REGIONS row.
        If values are empty, stores as "N/A".
        """
        try:
            parsed_data = {}
            
            # Parse CSV file
            lines = f['f'].splitlines()
            if len(lines) < 2:
                log.warning(f"Hi-C distance vs contact file {f.get('fn', 'unknown')} has insufficient lines")
                return
            
            # Parse header
            header = [col.strip() for col in lines[0].split(',')]
            if 'slope' not in header or 'mse' not in header:
                log.warning(f"Hi-C distance vs contact file {f.get('fn', 'unknown')} missing 'slope' or 'mse' column")
                return
            
            slope_idx = header.index('slope')
            mse_idx = header.index('mse')
            region_idx = header.index('region') if 'region' in header else None
            
            # Find ALL_REGIONS row (or first data row if no region column)
            found_row = False
            for line_idx, line in enumerate(lines[1:], start=1):
                if not line.strip():
                    continue
                
                data_row = [val.strip() for val in line.split(',')]
                if len(data_row) <= max(slope_idx, mse_idx):
                    continue
                
                # Check if this is the ALL_REGIONS row (if region column exists)
                if region_idx is not None:
                    if region_idx >= len(data_row) or data_row[region_idx].upper() != 'ALL_REGIONS':
                        continue
                
                # Extract slope
                slope_val = data_row[slope_idx] if slope_idx < len(data_row) else ''
                if slope_val and slope_val.strip():
                    try:
                        parsed_data['slope'] = float(slope_val)
                    except ValueError:
                        parsed_data['slope'] = "N/A"
                else:
                    parsed_data['slope'] = "N/A"
                
                # Extract mse
                mse_val = data_row[mse_idx] if mse_idx < len(data_row) else ''
                if mse_val and mse_val.strip():
                    try:
                        parsed_data['mse'] = float(mse_val)
                    except ValueError:
                        parsed_data['mse'] = "N/A"
                else:
                    parsed_data['mse'] = "N/A"
                
                # Found the row
                found_row = True
                break
            
            if not found_row:
                log.warning(f"Hi-C distance vs contact file {f.get('fn', 'unknown')}: Could not find ALL_REGIONS row")
                return
            
            if parsed_data:
                self.hic_dist_contact_data[f['s_name']] = parsed_data
                log.debug(f"Parsed Hi-C distance vs contact for '{f['s_name']}': {parsed_data}")
            else:
                log.warning(f"Hi-C distance vs contact file {f.get('fn', 'unknown')} has no valid data")
                
        except Exception as e:
            log.warning(f"Error parsing Hi-C distance vs contact file {f.get('fn', 'unknown')}: {e}")

    def parse_hic_loop_counts_tsv(self, f):
        """ Parses the Hi-C loop counts TSV file (loop_counts.tsv)
        
        Expected format:
        resolution      loop_count
        2000    0
        5000    0
        ...
        total   0
        
        Extracts the loop_count value from the row where resolution is "total".
        """
        try:
            parsed_data = {}
            
            # Parse TSV file
            lines = f['f'].splitlines()
            if len(lines) < 2:
                log.warning(f"Hi-C loop counts file {f.get('fn', 'unknown')} has insufficient lines")
                return
            
            # Parse header
            header = lines[0].strip().split('\t')
            if len(header) < 2:
                # Try splitting by multiple spaces
                header = lines[0].strip().split()
            
            if 'resolution' not in header or 'loop_count' not in header:
                log.warning(f"Hi-C loop counts file {f.get('fn', 'unknown')} missing 'resolution' or 'loop_count' column")
                return
            
            resolution_idx = header.index('resolution')
            loop_count_idx = header.index('loop_count')
            
            # Find the row where resolution is "total"
            found_total = False
            for line_idx, line in enumerate(lines[1:], start=1):
                if not line.strip():
                    continue
                
                # Split by tab first, then by whitespace if needed
                data_row = line.strip().split('\t')
                if len(data_row) < 2:
                    data_row = line.strip().split()
                
                if len(data_row) <= max(resolution_idx, loop_count_idx):
                    continue
                
                # Check if this is the "total" row
                resolution_val = data_row[resolution_idx].strip() if resolution_idx < len(data_row) else ''
                if resolution_val.lower() == 'total':
                    # Extract loop_count
                    loop_count_val = data_row[loop_count_idx].strip() if loop_count_idx < len(data_row) else ''
                    if loop_count_val:
                        try:
                            parsed_data['total_loops'] = int(loop_count_val)
                            found_total = True
                            log.debug(f"Found total loop count for '{f['s_name']}': {parsed_data['total_loops']}")
                        except ValueError:
                            log.warning(f"Hi-C loop counts file {f.get('fn', 'unknown')}: Invalid loop_count value '{loop_count_val}'")
                    break
            
            if not found_total:
                log.warning(f"Hi-C loop counts file {f.get('fn', 'unknown')}: Could not find 'total' row")
                return
            
            if parsed_data:
                self.hic_loop_counts_data[f['s_name']] = parsed_data
                log.debug(f"Parsed Hi-C loop counts for '{f['s_name']}': {parsed_data}")
            else:
                log.warning(f"Hi-C loop counts file {f.get('fn', 'unknown')} has no valid data")
                
        except Exception as e:
            log.warning(f"Error parsing Hi-C loop counts file {f.get('fn', 'unknown')}: {e}")

    def parse_hic_library_complexity_tsv(self, f):
        """ Parses the Hi-C library complexity TSV file
        
        Expected format:
        library N_total C       N_optical_rate  N_optical       N_pcr   U_observed
        test1-1 13801   59877.98        0.106876        1475    12326   11140.09
        
        Extracts the C column value (library complexity) for each library.
        """
        try:
            parsed_data = {}
            
            # Parse TSV file
            lines = f['f'].splitlines()
            if len(lines) < 2:
                log.warning(f"Hi-C library complexity file {f.get('fn', 'unknown')} has insufficient lines")
                return
            
            # Parse header
            header = lines[0].strip().split('\t')
            if len(header) < 2:
                # Try splitting by multiple spaces
                header = lines[0].strip().split()
            
            if 'library' not in header or 'C' not in header:
                log.warning(f"Hi-C library complexity file {f.get('fn', 'unknown')} missing 'library' or 'C' column")
                return
            
            library_idx = header.index('library')
            c_idx = header.index('C')
            
            # Parse data rows - try to match library name with sample name
            sample_name = f['s_name']
            found_match = False
            
            for line_idx, line in enumerate(lines[1:], start=1):
                if not line.strip():
                    continue
                
                # Split by tab first, then by whitespace if needed
                data_row = line.strip().split('\t')
                if len(data_row) < 2:
                    data_row = line.strip().split()
                
                if len(data_row) <= max(library_idx, c_idx):
                    continue
                
                # Extract library name and C value
                library_name = data_row[library_idx].strip() if library_idx < len(data_row) else ''
                c_val = data_row[c_idx].strip() if c_idx < len(data_row) else ''
                
                if library_name and c_val:
                    # Try to match library name with sample name (handle variations)
                    # Check if library name matches sample name (exact or partial match)
                    if library_name == sample_name or sample_name.endswith(library_name) or library_name in sample_name:
                        try:
                            complexity_value = float(c_val)
                            parsed_data['library_complexity'] = complexity_value
                            found_match = True
                            log.debug(f"Found library complexity for '{sample_name}' (matched library '{library_name}'): {complexity_value}")
                            break
                        except ValueError:
                            log.warning(f"Hi-C library complexity file {f.get('fn', 'unknown')}: Invalid C value '{c_val}'")
            
            # If no match found, use the first row (in case file is already split by sample)
            if not found_match and len(lines) > 1:
                first_data_line = lines[1].strip()
                if first_data_line:
                    data_row = first_data_line.split('\t')
                    if len(data_row) < 2:
                        data_row = first_data_line.split()
                    if len(data_row) > c_idx:
                        c_val = data_row[c_idx].strip() if c_idx < len(data_row) else ''
                        if c_val:
                            try:
                                complexity_value = float(c_val)
                                parsed_data['library_complexity'] = complexity_value
                                log.debug(f"Using first row library complexity for '{sample_name}': {complexity_value}")
                            except ValueError:
                                pass
            
            if parsed_data:
                self.hic_library_complexity_data[f['s_name']] = parsed_data
                log.debug(f"Parsed Hi-C library complexity for '{f['s_name']}': {parsed_data}")
            else:
                log.warning(f"Hi-C library complexity file {f.get('fn', 'unknown')} has no valid data")
                
        except Exception as e:
            log.warning(f"Error parsing Hi-C library complexity file {f.get('fn', 'unknown')}: {e}")

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
        
        # 14. Define Headers for Hi-C BAM Mapping Stats
        hic_mapstat_headers = OrderedDict()
        hic_mapstat_headers['mapped_rate'] = {
            'title': 'Hi-C Mapped Rate',
            'description': 'Hi-C: Mapped reads rate (mapped_1/total_1)',
            'min': 0,
            'max': 1,
            'format': '{:.4f}',
            'scale': 'RdYlGn'
        }
        hic_mapstat_headers['mapped_1'] = {
            'title': 'Hi-C Mapped',
            'description': 'Hi-C: Number of mapped reads',
            'format': '{:,.0f}',
            'scale': 'Blues'
        }
        hic_mapstat_headers['total_1'] = {
            'title': 'Hi-C Total',
            'description': 'Hi-C: Total number of reads',
            'format': '{:,.0f}',
            'scale': 'Blues'
        }
        
        # 15. Define Headers for Hi-C Pair Alignment Stats
        hic_pairstat_headers = OrderedDict()
        hic_pairstat_headers['Unique_paired_alignments'] = {
            'title': 'Hi-C Unique Pairs',
            'description': 'Hi-C: Number of unique paired alignments',
            'format': '{:,.0f}',
            'scale': 'Greens'
        }
        hic_pairstat_headers['Unique_paired_alignments_percent'] = {
            'title': 'Hi-C Unique %',
            'description': 'Hi-C: Percentage of unique paired alignments',
            'min': 0,
            'max': 100,
            'suffix': '%',
            'format': '{:.2f}',
            'scale': 'Greens'
        }
        hic_pairstat_headers['Total_pairs_processed'] = {
            'title': 'Hi-C Total Pairs',
            'description': 'Hi-C: Total pairs processed',
            'format': '{:,.0f}',
            'scale': 'Blues'
        }
        hic_pairstat_headers['Mapped_ratio'] = {
            'title': 'Hi-C Mapped Ratio',
            'description': 'Hi-C: Mapped ratio (1 - Unmapped_pairs/Total_pairs_processed)',
            'min': 0,
            'max': 1,
            'format': '{:.4f}',
            'scale': 'RdYlGn'
        }
        hic_pairstat_headers['unique_ratio'] = {
            'title': 'Hi-C Unique Ratio',
            'description': 'Hi-C: Unique paired alignments ratio (Unique_paired_alignments/Total_pairs_processed)',
            'min': 0,
            'max': 1,
            'format': '{:.4f}',
            'scale': 'RdYlGn'
        }
        
        # 16. Define Headers for Hi-C Processing Stats (RSstat)
        hic_rsstat_headers = OrderedDict()
        hic_rsstat_headers['Valid_interaction_pairs'] = {
            'title': 'Hi-C Valid Pairs',
            'description': 'Hi-C: Number of valid interaction pairs',
            'format': '{:,.0f}',
            'scale': 'Greens'
        }
        hic_rsstat_headers['Dangling_end_pairs'] = {
            'title': 'Hi-C Dangling',
            'description': 'Hi-C: Number of dangling end pairs',
            'format': '{:,.0f}',
            'scale': 'Oranges'
        }
        hic_rsstat_headers['Religation_pairs'] = {
            'title': 'Hi-C Religation',
            'description': 'Hi-C: Number of religation pairs',
            'format': '{:,.0f}',
            'scale': 'Oranges'
        }
        hic_rsstat_headers['dangling_rate'] = {
            'title': 'Hi-C Dangling Rate',
            'description': 'Hi-C: Dangling end pairs rate (Dangling_end_pairs/Total_pairs_processed)',
            'min': 0,
            'max': 1,
            'format': '{:.4f}',
            'scale': 'OrRd'
        }
        hic_rsstat_headers['self_ligation_rate'] = {
            'title': 'Hi-C Self-Ligation Rate',
            'description': 'Hi-C: Self-ligation rate ((Self_Cycle_pairs+Religation_pairs)/Total_pairs_processed)',
            'min': 0,
            'max': 1,
            'format': '{:.4f}',
            'scale': 'OrRd'
        }
        
        # 17. Define Headers for Hi-C Deduplication Stats
        hic_dedup_headers = OrderedDict()
        hic_dedup_headers['frac_dups'] = {
            'title': 'Hi-C Frac Dups',
            'description': 'Hi-C: Fraction of duplicates (from summary/frac_dups)',
            'min': 0,
            'max': 1,
            'format': '{:.4f}',
            'scale': 'OrRd'
        }
        hic_dedup_headers['frac_cis2k'] = {
            'title': 'Hi-C Frac Cis 2kb',
            'description': 'Hi-C: Fraction of cis interactions within 2kb (cis_2kb/Total_pairs_processed)',
            'min': 0,
            'max': 1,
            'format': '{:.4f}',
            'scale': 'Blues'
        }
        hic_dedup_headers['frac_trans'] = {
            'title': 'Hi-C Frac Trans',
            'description': 'Hi-C: Fraction of trans interactions (trans/total)',
            'min': 0,
            'max': 1,
            'format': '{:.4f}',
            'scale': 'Blues'
        }
        hic_dedup_headers['cis_trans_ratio'] = {
            'title': 'Hi-C Cis/Trans Ratio',
            'description': 'Hi-C: Ratio of cis to trans interactions (cis/trans)',
            'min': 0,
            'format': '{:.4f}',
            'scale': 'RdYlGn'
        }
        
        # 18. Define Headers for Hi-C Distance vs Contact Stats
        hic_dist_contact_headers = OrderedDict()
        hic_dist_contact_headers['slope'] = {
            'title': 'Hi-C Slope',
            'description': 'Hi-C: Distance vs contact slope (log-log fit)',
            'format': '{:.4f}',
            'scale': 'RdYlGn'
        }
        hic_dist_contact_headers['mse'] = {
            'title': 'Hi-C MSE',
            'description': 'Hi-C: Mean squared error of distance vs contact fit',
            'format': '{:.4f}',
            'scale': 'OrRd',
            'min': 0
        }
        
        # 19. Define Headers for Hi-C Loop Counts Stats
        hic_loop_counts_headers = OrderedDict()
        hic_loop_counts_headers['total_loops'] = {
            'title': 'Hi-C Total Loops',
            'description': 'Hi-C: Total number of loops detected across all resolutions',
            'format': '{:,.0f}',
            'scale': 'Blues',
            'min': 0
        }
        
        # 20. Define Headers for Hi-C Library Complexity Stats
        hic_library_complexity_headers = OrderedDict()
        hic_library_complexity_headers['library_complexity'] = {
            'title': 'Hi-C Library Complexity',
            'description': 'Hi-C: Library complexity (C)',
            'format': '{:,.2f}',
            'scale': 'RdYlGn',
            'min': 0
        }

        # Calculate cross-file metrics before adding to general stats
        # Match data across files by sample name
        all_samples = set()
        all_samples.update(self.hic_pairstat_data.keys())
        all_samples.update(self.hic_rsstat_data.keys())
        all_samples.update(self.hic_dedup_stats_data.keys())
        
        for sample_name in all_samples:
            # Get Total_pairs_processed from pairstat data
            total_pairs_processed = None
            if sample_name in self.hic_pairstat_data:
                total_pairs_processed = self.hic_pairstat_data[sample_name].get('Total_pairs_processed')
            
            # Calculate RSstat metrics if Total_pairs_processed is available
            if sample_name in self.hic_rsstat_data and total_pairs_processed is not None and total_pairs_processed > 0:
                rsstat_data = self.hic_rsstat_data[sample_name]
                
                # Calculate dangling_rate = Dangling_end_pairs / Total_pairs_processed
                if 'Dangling_end_pairs' in rsstat_data:
                    dangling_pairs = rsstat_data['Dangling_end_pairs']
                    dangling_rate = float(dangling_pairs) / float(total_pairs_processed)
                    rsstat_data['dangling_rate'] = dangling_rate
                    log.debug(f"Calculated dangling_rate for '{sample_name}': {dangling_rate:.4f}")
                
                # Calculate self_ligation_rate = (Self_Cycle_pairs + Religation_pairs) / Total_pairs_processed
                self_cycle = rsstat_data.get('Self_Cycle_pairs', 0)
                religation = rsstat_data.get('Religation_pairs', 0)
                self_ligation_rate = float(self_cycle + religation) / float(total_pairs_processed)
                rsstat_data['self_ligation_rate'] = self_ligation_rate
                log.debug(f"Calculated self_ligation_rate for '{sample_name}': {self_ligation_rate:.4f}")
            
            # Calculate dedup metrics if Total_pairs_processed is available
            if sample_name in self.hic_dedup_stats_data and total_pairs_processed is not None and total_pairs_processed > 0:
                dedup_data = self.hic_dedup_stats_data[sample_name]
                
                # Calculate frac_cis2k = cis_2kb / Total_pairs_processed
                cis_2kb = dedup_data.get('cis_2kb')
                if cis_2kb is not None:
                    frac_cis2k = float(cis_2kb) / float(total_pairs_processed)
                    dedup_data['frac_cis2k'] = frac_cis2k
                    log.debug(f"Calculated frac_cis2k for '{sample_name}': {frac_cis2k:.4f}")

        # 17. Add to General Stats
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
        self.general_stats_addcols(self.hic_mapstat_data, hic_mapstat_headers)
        self.general_stats_addcols(self.hic_pairstat_data, hic_pairstat_headers)
        self.general_stats_addcols(self.hic_rsstat_data, hic_rsstat_headers)
        self.general_stats_addcols(self.hic_dedup_stats_data, hic_dedup_headers)
        self.general_stats_addcols(self.hic_dist_contact_data, hic_dist_contact_headers)
        self.general_stats_addcols(self.hic_loop_counts_data, hic_loop_counts_headers)
        self.general_stats_addcols(self.hic_library_complexity_data, hic_library_complexity_headers)

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