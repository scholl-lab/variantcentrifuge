# Output and Reporting Stages Analysis

## Current State Analysis

### Code Location
- **File**: `pipeline.py`
- **Lines**:
  - 2342-2399: Variant ID generation
  - 2446-2481: Named column addition
  - 2484-2510: Final filtering
  - 2511-2562: Pseudonymization
  - 2563-2592: Final output writing
  - 2649-2773: Excel generation
  - 2736-2768: HTML report generation
  - 2679-2734: IGV report generation
- **External Files**:
  - `converter.py`: TSV to Excel conversion
  - `generate_html_report.py`: HTML report generation
  - `generate_igv_report.py`: IGV integration
  - `pseudonymizer.py`: Sample pseudonymization

### Current Functionality

The output phase handles multiple output formats and post-processing:

#### 1. Post-Processing
- Add variant identifiers (VAR_ID)
- Add named blank columns
- Apply final filtering
- Pseudonymize sample names

#### 2. Output Formats
- Primary TSV output
- Excel workbook with multiple sheets
- Interactive HTML report
- IGV.js integration files
- Metadata files

#### 3. Report Generation
- Summary statistics
- Interactive tables
- Genomic visualizations
- External database links

### Current Issues

1. **Sequential Processing**: Each output format generated separately
2. **File I/O Overhead**: Multiple reads of the same data
3. **Memory Usage**: Some operations load entire file
4. **Complex Dependencies**: Output generation order matters
5. **Error Recovery**: Difficult to regenerate single output

## Refactored Design

### Core Output Stages

```python
# stages/output.py

class VariantIdentifierStage(Stage):
    """Add unique variant identifiers."""
    
    @property
    def name(self) -> str:
        return "variant_identifier"
    
    @property
    def dependencies(self) -> Set[str]:
        # Run after all analysis
        return {"statistics_generation"} or {"chunked_analysis"}
    
    def _process(self, context: PipelineContext) -> PipelineContext:
        """Add VAR_ID column."""
        if context.use_chunked_processing:
            # Process file with streaming
            self._add_identifiers_streaming(context)
        else:
            # Process DataFrame
            df = self._add_identifiers_dataframe(context.current_dataframe)
            context.current_dataframe = df
        
        return context
    
    def _add_identifiers_dataframe(self, df: pd.DataFrame) -> pd.DataFrame:
        """Add VAR_ID to DataFrame."""
        import hashlib
        
        # Generate IDs
        var_ids = []
        for idx, row in df.iterrows():
            variant_str = f"{row['CHROM']}{row['POS']}{row['REF']}{row['ALT']}"
            short_hash = hashlib.md5(variant_str.encode()).hexdigest()[:4]
            var_id = f"var_{idx+1:04d}_{short_hash}"
            var_ids.append(var_id)
        
        # Insert at beginning
        df.insert(0, 'VAR_ID', var_ids)
        return df
    
    def _add_identifiers_streaming(self, context: PipelineContext) -> None:
        """Add VAR_ID with streaming for large files."""
        input_file = context.current_file
        output_file = context.intermediate_dir / f"{context.base_name}.with_ids.tsv.gz"
        
        with smart_open(input_file, 'r') as inp, \
             smart_open(output_file, 'w') as out:
            
            # Process header
            header = inp.readline().strip()
            header_parts = header.split('\t')
            
            # Find required columns
            col_indices = {}
            for col in ['CHROM', 'POS', 'REF', 'ALT']:
                if col in header_parts:
                    col_indices[col] = header_parts.index(col)
            
            # Write new header
            out.write('VAR_ID\t' + header + '\n')
            
            # Process data
            for idx, line in enumerate(inp, start=1):
                fields = line.strip().split('\t')
                
                # Generate ID
                variant_parts = []
                for col in ['CHROM', 'POS', 'REF', 'ALT']:
                    if col in col_indices:
                        variant_parts.append(fields[col_indices[col]])
                
                variant_str = ''.join(variant_parts)
                short_hash = hashlib.md5(variant_str.encode()).hexdigest()[:4]
                var_id = f"var_{idx:04d}_{short_hash}"
                
                out.write(f"{var_id}\t{line}")
        
        context.current_file = output_file


class FinalFilteringStage(Stage):
    """Apply final pandas-style filtering."""
    
    @property
    def name(self) -> str:
        return "final_filtering"
    
    @property
    def dependencies(self) -> Set[str]:
        return {"variant_identifier"}
    
    def _process(self, context: PipelineContext) -> PipelineContext:
        """Apply final filter if configured."""
        final_filter = context.config.get('final_filter')
        if not final_filter:
            return context
        
        if context.use_chunked_processing:
            # Convert to DataFrame for filtering
            df = pd.read_csv(
                context.current_file,
                sep='\t',
                dtype=str,
                keep_default_na=False
            )
        else:
            df = context.current_dataframe
        
        # Apply filter
        from ..filters import filter_dataframe_with_query
        filtered_df = filter_dataframe_with_query(df, final_filter)
        
        if len(filtered_df) == 0:
            logger.warning("No variants passed final filter")
        
        context.current_dataframe = filtered_df
        context.use_chunked_processing = False  # Now in memory
        
        return context


class PseudonymizationStage(Stage):
    """Pseudonymize sample identifiers."""
    
    @property
    def name(self) -> str:
        return "pseudonymization"
    
    @property
    def dependencies(self) -> Set[str]:
        return {"final_filtering"}
    
    def _process(self, context: PipelineContext) -> PipelineContext:
        """Apply pseudonymization if configured."""
        if not context.config.get('pseudonymize'):
            return context
        
        # Extract sample list from data
        samples = self._extract_samples(context)
        
        # Get pedigree data if available
        ped_data = context.get_stage_result('pedigree_data')
        
        # Apply pseudonymization
        from ..pseudonymizer import apply_pseudonymization
        
        if context.current_dataframe is not None:
            # DataFrame mode
            df_lines = self._dataframe_to_lines(context.current_dataframe)
            result_lines, pseudonymizer = apply_pseudonymization(
                df_lines, samples, context.config, ped_data
            )
            context.current_dataframe = self._lines_to_dataframe(result_lines)
        else:
            # File mode
            raise NotImplementedError("Streaming pseudonymization not yet implemented")
        
        # Save mapping
        context.pseudonymizer = pseudonymizer
        self._save_pseudonymization_files(context, pseudonymizer)
        
        return context
    
    def _save_pseudonymization_files(self, context: PipelineContext, 
                                    pseudonymizer: Any) -> None:
        """Save pseudonymization mapping and PED if requested."""
        if context.config.get('pseudonymize_table'):
            # Save in parent directory of output
            output_parent = context.output_dir.parent
            table_path = output_parent / Path(context.config['pseudonymize_table']).name
            pseudonymizer.save_mapping(table_path)
            logger.info(f"Pseudonymization mapping saved to: {table_path}")
            
            # Create pseudonymized PED if requested
            if context.config.get('pseudonymize_ped') and context.config.get('ped_file'):
                ped_output = table_path.with_suffix('.ped')
                pseudonymizer.pseudonymize_ped_file(
                    context.config['ped_file'], 
                    ped_output
                )
                logger.info(f"Pseudonymized PED saved to: {ped_output}")


class TSVOutputStage(Stage):
    """Write final TSV output."""
    
    @property
    def name(self) -> str:
        return "tsv_output"
    
    @property
    def dependencies(self) -> Set[str]:
        return {"pseudonymization"}
    
    def _process(self, context: PipelineContext) -> PipelineContext:
        """Write final TSV file."""
        # Determine output path
        if context.args.output_file in ['stdout', '-']:
            # Write to stdout
            self._write_to_stdout(context)
            context.final_output_path = None
        else:
            # Write to file
            output_path = context.output_dir / Path(context.args.output_file).name
            self._write_to_file(context, output_path)
            context.final_output_path = output_path
        
        return context
    
    def _write_to_stdout(self, context: PipelineContext) -> None:
        """Write output to stdout."""
        if context.current_dataframe is not None:
            context.current_dataframe.to_csv(
                sys.stdout,
                sep='\t',
                index=False,
                na_rep=''
            )
        else:
            # Stream from file
            with smart_open(context.current_file, 'r') as f:
                for line in f:
                    sys.stdout.write(line)
    
    def _write_to_file(self, context: PipelineContext, output_path: Path) -> None:
        """Write output to file."""
        if context.current_dataframe is not None:
            context.current_dataframe.to_csv(
                output_path,
                sep='\t',
                index=False,
                na_rep=''
            )
        else:
            # Copy from current file
            shutil.copy2(context.current_file, output_path)
        
        logger.info(f"Output saved to: {output_path}")


class ExcelReportStage(Stage):
    """Generate Excel workbook with multiple sheets."""
    
    @property
    def name(self) -> str:
        return "excel_report"
    
    @property
    def dependencies(self) -> Set[str]:
        return {"tsv_output", "metadata_generation"}
    
    @property
    def can_run_parallel(self) -> bool:
        return True  # Independent of other reports
    
    def _process(self, context: PipelineContext) -> PipelineContext:
        """Generate Excel report if requested."""
        if not context.args.xlsx or not context.final_output_path:
            return context
        
        excel_path = context.final_output_path.with_suffix('.xlsx')
        
        # Create workbook
        from ..converter import ExcelReportGenerator
        
        generator = ExcelReportGenerator(excel_path)
        
        # Add main results
        generator.add_tsv_sheet(
            context.final_output_path,
            sheet_name='Variants'
        )
        
        # Add metadata
        if context.metadata_path.exists():
            generator.add_tsv_sheet(
                context.metadata_path,
                sheet_name='Metadata'
            )
        
        # Add statistics
        stats_file = context.config.get('stats_output_file')
        if stats_file and Path(stats_file).exists():
            generator.add_tsv_sheet(
                Path(stats_file),
                sheet_name='Statistics'
            )
        
        # Add gene burden if available
        gene_burden_file = context.output_dir / f"{context.base_name}.gene_burden.tsv"
        if gene_burden_file.exists():
            generator.add_tsv_sheet(
                gene_burden_file,
                sheet_name='Gene Burden'
            )
        
        # Finalize with formatting
        generator.finalize(context.config)
        
        context.excel_path = excel_path
        logger.info(f"Excel report saved to: {excel_path}")
        
        return context


class HTMLReportStage(Stage):
    """Generate interactive HTML report."""
    
    @property
    def name(self) -> str:
        return "html_report"
    
    @property
    def dependencies(self) -> Set[str]:
        return {"tsv_output", "igv_report"}  # Needs IGV links if enabled
    
    @property
    def can_run_parallel(self) -> bool:
        return True
    
    def _process(self, context: PipelineContext) -> PipelineContext:
        """Generate HTML report if requested."""
        if not context.args.html_report or not context.final_output_path:
            return context
        
        report_dir = context.output_dir / 'report'
        report_dir.mkdir(exist_ok=True)
        
        # First generate JSON data
        from ..converter import produce_report_json
        produce_report_json(context.final_output_path, context.output_dir)
        
        # Then generate HTML
        from ..generate_html_report import generate_html_report
        
        variants_json = report_dir / 'variants.json'
        summary_json = report_dir / 'summary.json'
        
        generate_html_report(
            variants_json=variants_json,
            summary_json=summary_json,
            output_dir=report_dir,
            cfg=context.config
        )
        
        context.html_report_path = report_dir / 'index.html'
        logger.info(f"HTML report saved to: {context.html_report_path}")
        
        return context


class IGVReportStage(Stage):
    """Generate IGV.js visualization files."""
    
    @property
    def name(self) -> str:
        return "igv_report"
    
    @property
    def dependencies(self) -> Set[str]:
        return {"tsv_output"}
    
    @property
    def can_run_parallel(self) -> bool:
        return True
    
    def _process(self, context: PipelineContext) -> PipelineContext:
        """Generate IGV reports if enabled."""
        if not context.config.get('igv_enabled') or not context.final_output_path:
            return context
        
        # Validate IGV configuration
        self._validate_igv_config(context)
        
        report_dir = context.output_dir / 'report'
        report_dir.mkdir(exist_ok=True)
        
        from ..generate_igv_report import generate_igv_report
        
        generate_igv_report(
            variants_tsv=context.final_output_path,
            output_dir=report_dir,
            bam_mapping_file=context.config.get('bam_mapping_file'),
            igv_reference=context.config.get('igv_reference'),
            integrate_into_main=True,
            igv_fasta=context.config.get('igv_fasta'),
            igv_ideogram=context.config.get('igv_ideogram'),
            igv_flanking=context.config.get('igv_flanking', 50)
        )
        
        context.igv_report_generated = True
        logger.info("IGV reports and mapping file generated")
        
        return context
    
    def _validate_igv_config(self, context: PipelineContext) -> None:
        """Validate IGV configuration."""
        if not context.config.get('bam_mapping_file'):
            raise ValueError("IGV integration requires --bam-mapping-file")
        
        if not context.config.get('igv_reference') and not context.config.get('igv_fasta'):
            raise ValueError("IGV integration requires either --igv-reference or --igv-fasta")


class MetadataGenerationStage(Stage):
    """Generate metadata file with run information."""
    
    @property
    def name(self) -> str:
        return "metadata_generation"
    
    @property
    def dependencies(self) -> Set[str]:
        return {"tsv_output"}
    
    @property
    def can_run_parallel(self) -> bool:
        return True
    
    def _process(self, context: PipelineContext) -> PipelineContext:
        """Generate metadata file."""
        metadata_path = context.output_dir / f"{context.base_name}.metadata.tsv"
        
        with open(metadata_path, 'w') as f:
            f.write("Parameter\tValue\n")
            
            # Basic information
            self._write_metadata(f, "Tool", "variantcentrifuge")
            self._write_metadata(f, "Version", context.config.get('version', 'N/A'))
            self._write_metadata(f, "Run_start_time", context.start_time.isoformat())
            self._write_metadata(f, "Run_end_time", datetime.now().isoformat())
            
            # Duration
            duration = (datetime.now() - context.start_time).total_seconds()
            self._write_metadata(f, "Run_duration_seconds", str(duration))
            
            # Command line
            self._write_metadata(f, "Command_line", ' '.join(sys.argv))
            
            # Configuration
            for key, value in context.config.items():
                if not key.startswith('_'):  # Skip internal keys
                    self._write_metadata(f, f"config.{key}", str(value))
            
            # Tool versions
            self._write_tool_versions(f)
        
        context.metadata_path = metadata_path
        logger.info(f"Metadata saved to: {metadata_path}")
        
        return context
    
    def _write_metadata(self, file, key: str, value: str) -> None:
        """Write sanitized metadata entry."""
        from ..utils import sanitize_metadata_field
        key = sanitize_metadata_field(key)
        value = sanitize_metadata_field(value)
        file.write(f"{key}\t{value}\n")
    
    def _write_tool_versions(self, file) -> None:
        """Write external tool versions."""
        from ..utils import get_tool_version
        
        tools = ['snpEff', 'bcftools', 'SnpSift', 'bedtools']
        for tool in tools:
            version = get_tool_version(tool)
            self._write_metadata(file, f"tool.{tool}_version", version)
```

### Parallel Report Generation

```python
# stages/output_parallel.py

class ParallelReportGenerationStage(Stage):
    """Generate all reports in parallel."""
    
    @property
    def name(self) -> str:
        return "parallel_report_generation"
    
    @property
    def dependencies(self) -> Set[str]:
        return {"tsv_output"}
    
    def _process(self, context: PipelineContext) -> PipelineContext:
        """Generate all requested reports in parallel."""
        report_stages = []
        
        # Add requested report stages
        if context.args.xlsx:
            report_stages.append(ExcelReportStage())
        
        if context.config.get('igv_enabled'):
            report_stages.append(IGVReportStage())
        
        if context.args.html_report:
            # HTML depends on IGV, so handle specially
            if context.config.get('igv_enabled'):
                # Will run after IGV
                pass
            else:
                report_stages.append(HTMLReportStage())
        
        # Always generate metadata
        report_stages.append(MetadataGenerationStage())
        
        # Run in parallel
        with ThreadPoolExecutor(max_workers=len(report_stages)) as executor:
            futures = {
                executor.submit(stage, context): stage 
                for stage in report_stages
            }
            
            for future in as_completed(futures):
                stage = futures[future]
                try:
                    result = future.result()
                    logger.info(f"Completed: {stage.name}")
                except Exception as e:
                    logger.error(f"Failed to generate {stage.name}: {e}")
                    raise
        
        # Run HTML report after IGV if both requested
        if context.args.html_report and context.config.get('igv_enabled'):
            HTMLReportStage()(context)
        
        return context
```

## Testing Strategy

### Unit Tests

```python
# tests/unit/stages/test_output.py

class TestOutputStages:
    
    def test_variant_identifier_generation(self, sample_df):
        """Test VAR_ID generation."""
        context = create_test_context()
        context.current_dataframe = sample_df
        
        stage = VariantIdentifierStage()
        result = stage(context)
        
        df = result.current_dataframe
        assert 'VAR_ID' in df.columns
        assert df.columns[0] == 'VAR_ID'  # First column
        
        # Check format
        assert all(df['VAR_ID'].str.match(r'^var_\d{4}_[a-f0-9]{4}$'))
        
        # Check uniqueness
        assert len(df['VAR_ID'].unique()) == len(df)
    
    def test_pseudonymization(self, sample_df_with_samples):
        """Test sample name pseudonymization."""
        context = create_test_context()
        context.current_dataframe = sample_df_with_samples
        context.config['pseudonymize'] = True
        context.config['pseudonymize_schema'] = 'sequential'
        
        stage = PseudonymizationStage()
        result = stage(context)
        
        df = result.current_dataframe
        
        # Check samples replaced
        assert 'SAMPLE_001' in df['GT'].iloc[0]
        assert 'Sample1' not in df['GT'].iloc[0]
        
        # Check mapping saved
        assert context.pseudonymizer is not None
    
    def test_excel_generation(self, tmp_path):
        """Test Excel report generation."""
        # Create test data
        tsv_file = tmp_path / "test.tsv"
        tsv_file.write_text("CHROM\tPOS\tREF\tALT\nchr1\t100\tA\tG\n")
        
        context = create_test_context()
        context.final_output_path = tsv_file
        context.args.xlsx = True
        
        stage = ExcelReportStage()
        result = stage(context)
        
        assert result.excel_path.exists()
        assert result.excel_path.suffix == '.xlsx'
        
        # Verify content
        df = pd.read_excel(result.excel_path, sheet_name='Variants')
        assert len(df) == 1
        assert df.iloc[0]['CHROM'] == 'chr1'
```

### Integration Tests

```python
# tests/integration/test_output_pipeline.py

def test_complete_output_generation(test_analysis_results):
    """Test all output formats generated correctly."""
    context = create_test_context(
        args=create_test_args(
            xlsx=True,
            html_report=True,
            output_file="results.tsv"
        ),
        config={
            'igv_enabled': True,
            'bam_mapping_file': 'bams.txt',
            'igv_reference': 'hg38'
        }
    )
    
    # Set up test data
    context.current_dataframe = test_analysis_results
    
    # Run output stages
    stages = [
        VariantIdentifierStage(),
        FinalFilteringStage(),
        PseudonymizationStage(),
        TSVOutputStage(),
        ParallelReportGenerationStage()
    ]
    
    runner = PipelineRunner()
    result = runner.run(stages, context)
    
    # Verify all outputs
    assert result.final_output_path.exists()
    assert result.excel_path.exists()
    assert result.html_report_path.exists()
    assert result.metadata_path.exists()
    assert result.igv_report_generated
```

## Performance Considerations

### Parallel Report Generation

Generate independent reports concurrently:
- Excel generation
- IGV report generation  
- Metadata generation
- HTML (if no IGV dependency)

### Streaming Operations

For large files:
- Variant ID addition
- Pseudonymization
- Final filtering (convert to DataFrame only when needed)

### Memory Optimization

- Don't reload data for each output format
- Share DataFrame between stages when possible
- Stream large files

## Migration Plan

### Week 1: Core Output Stages
1. Implement VariantIdentifierStage
2. Implement TSVOutputStage
3. Test basic output generation

### Week 2: Post-Processing
1. Implement FinalFilteringStage
2. Implement PseudonymizationStage
3. Test transformations

### Week 3: Report Generation
1. Implement ExcelReportStage
2. Implement HTMLReportStage
3. Implement IGVReportStage

### Week 4: Optimization
1. Implement ParallelReportGenerationStage
2. Add streaming support for large files
3. Performance testing

## Success Metrics

- All output formats identical to current implementation
- Report generation 2-3x faster with parallel execution
- Memory usage constant regardless of file size
- Each output format independently regeneratable
- Clear error messages for missing dependencies