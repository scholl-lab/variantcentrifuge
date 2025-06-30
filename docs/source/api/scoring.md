# Scoring Module

The `scoring` module provides functionality for applying custom variant scoring formulas to variant data.

## Functions

### read_scoring_config

```python
def read_scoring_config(config_path: str) -> Dict[str, Any]
```

Read and parse the scoring configuration files from a directory.

**Parameters:**
- `config_path` (str): Path to the scoring configuration directory

**Returns:**
- `Dict[str, Any]`: A dictionary containing the parsed variables and formulas

**Raises:**
- `FileNotFoundError`: If configuration files are not found
- `json.JSONDecodeError`: If JSON parsing fails

**Expected Files:**
- `variable_assignment_config.json`: Maps DataFrame columns to formula variable names
- `formula_config.json`: Contains the scoring formulas

**Example:**
```python
from variantcentrifuge.scoring import read_scoring_config

config = read_scoring_config("scoring/nephro_variant_score")
print(config['variables'])  # Variable mappings
print(config['formulas'])   # Scoring formulas
```

### convert_to_numeric

```python
def convert_to_numeric(series: pd.Series, default: float = 0.0) -> pd.Series
```

Convert a pandas Series to numeric, handling empty strings and other non-numeric values.

**Parameters:**
- `series` (pd.Series): The series to convert
- `default` (float): The default value to use for non-numeric entries (default: 0.0)

**Returns:**
- `pd.Series`: The numeric series

**Example:**
```python
import pandas as pd
from variantcentrifuge.scoring import convert_to_numeric

# Handle mixed data types
data = pd.Series(['1.5', '', '2.0', '.', '3.5'])
numeric_data = convert_to_numeric(data, default=0.0)
# Result: [1.5, 0.0, 2.0, 0.0, 3.5]
```

### apply_scoring

```python
def apply_scoring(df: pd.DataFrame, scoring_config: Dict[str, Any]) -> pd.DataFrame
```

Apply scoring formulas to a DataFrame of variants.

**Parameters:**
- `df` (pd.DataFrame): The DataFrame containing annotated variant data
- `scoring_config` (Dict[str, Any]): The parsed scoring configuration

**Returns:**
- `pd.DataFrame`: The DataFrame with new score columns added

**Process:**
1. Maps original column names to formula variable names
2. Handles missing columns by creating them with default values
3. Converts numeric columns to proper numeric types
4. Evaluates formulas using pandas.eval()
5. Adds resulting scores as new columns

**Example:**
```python
import pandas as pd
from variantcentrifuge.scoring import read_scoring_config, apply_scoring

# Load variant data
variants_df = pd.read_csv("variants.tsv", sep="\t")

# Load scoring configuration
config = read_scoring_config("scoring/nephro_variant_score")

# Apply scoring
scored_df = apply_scoring(variants_df, config)

# Access the new score column
print(scored_df['nephro_variant_score'])
```

## Configuration Format

### Variable Assignment Configuration

The `variable_assignment_config.json` file maps VCF annotation fields to formula variables:

```json
{
  "variables": {
    "original_column_name": "formula_variable_name|default:value",
    "dbNSFP_gnomAD_exomes_AF": "gnomade_variant|default:0.0",
    "ANN[0].IMPACT": "impact_variant|default:''"
  }
}
```

### Formula Configuration

The `formula_config.json` file contains scoring formulas:

```json
{
  "formulas": [
    {
      "score_name": "pandas_eval_expression"
    }
  ]
}
```

## Formula Syntax

Formulas use pandas eval syntax and support:

- **Arithmetic operations**: `+`, `-`, `*`, `/`, `**`
- **Comparisons**: `==`, `!=`, `<`, `>`, `<=`, `>=`
- **Logical operations**: `&` (and), `|` (or), `~` (not)
- **Conditional logic**: `(condition) * value`
- **String operations**: `.str.contains()`, `.str.lower()`
- **Mathematical functions**: Available through numeric operations

**Example Formula:**
```python
"1 / (1 + 2.718281828459045 ** (-((intercept) + (var1 * coef1) + (var2 * coef2))))"
```

## Integration Example

```python
from variantcentrifuge.pipeline import run_pipeline
from variantcentrifuge.scoring import read_scoring_config

# Configuration with scoring
config = {
    "gene_name": "BRCA1",
    "vcf_file": "input.vcf.gz",
    "scoring_config_path": "scoring/nephro_variant_score",
    "output_file": "scored_variants.tsv"
}

# The pipeline will automatically:
# 1. Load the scoring configuration
# 2. Apply scoring after variant analysis
# 3. Include scores in the output
run_pipeline(config)
```

## Notes

- Missing columns are handled gracefully with default values
- Numeric columns are automatically converted from strings
- Column renaming persists in the output DataFrame
- Multiple formulas can be applied in a single configuration
- Formulas are evaluated using pandas' Python engine for maximum compatibility