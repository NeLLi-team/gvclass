# GVClass BLAST Improvement Plan

## Issue #1: Identity Percentage and Hit Sorting

### Problem Statement
The current BLAST implementation had two major issues:
1. Identity percentage was hardcoded to 100% instead of being calculated from actual alignment
2. BLAST hits were not properly sorted by score before selection for phylogenetic analysis

### Solution Implemented

#### 1. Identity Percentage Calculation
- Modified `run_blastp()` in `src/core/blast.py` to extract real identity from pyswrd's result object
- The pyswrd library provides alignment details through `hit.result` which includes:
  - `identity()`: Returns fraction of identical residues (0.0 to 1.0)
  - `query_start/end`: Alignment positions in query
  - `target_start/end`: Alignment positions in target
  - `alignment`: String representation of alignment

#### 2. Proper Hit Sorting by Score
- Collect all hits before writing to output file
- Sort hits by score (descending) and e-value (ascending) as secondary criterion
- This ensures best-scoring hits are prioritized for downstream analysis

#### 3. Enhanced parse_blastp Function
- Added configurable parameters:
  - `max_hits_per_query`: Limit hits per query sequence (default: 10)
  - `min_identity`: Minimum identity threshold (default: 30%)
- Maintains score-based ordering throughout processing
- Returns unique subject IDs sorted by score

### Technical Details

The BLAST output format now includes accurate metrics:
```
query_id  target_id  identity%  aln_len  mismatches  gaps  q_start  q_end  s_start  s_end  evalue  score
```

### Benefits
1. **More accurate phylogenetic placement**: Using actual identity percentages ensures only truly similar sequences are included
2. **Better reference selection**: High-scoring hits are prioritized, leading to more meaningful phylogenetic trees
3. **Configurable thresholds**: Users can adjust identity and hit count parameters based on their needs
4. **Improved performance**: Sorting and filtering reduce unnecessary sequences in alignment/tree building

### Files Modified
- `src/core/blast.py`: 
  - `run_blastp()`: Added identity calculation and hit sorting
  - `parse_blastp()`: Enhanced with score-based selection and filtering

### Testing Recommendations
1. Run with example data to verify identity percentages are realistic (not 100%)
2. Check that BLAST output is sorted by score (highest first)
3. Verify phylogenetic trees include appropriate reference sequences
4. Compare results with previous version to ensure improvements

### Future Enhancements
- Consider adding coverage percentage as additional filter
- Allow user-configurable e-value threshold (currently fixed at 1e-5)
- Add option to weight hits by both identity and score
- Implement adaptive thresholds based on query characteristics