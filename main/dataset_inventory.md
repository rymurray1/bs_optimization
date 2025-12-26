# Dataset Inventory - Composition and Process Data

## Purpose
Manual examination and cataloging of all datasets in composition_data.csv and process_data.csv

---

## COMPOSITION DATA BLOCKS

### Dataset 1: Ausenco-Tiger-Gold (Lines 3-8)
- **Project**: Ausenco-Tiger-Gold, December 10, 2025
- **Format**: Composite-level data
- **Columns**: Composite, Au (g/t), Ag (g/t), Fe (%), Cu (g/t), Pb (g/t), ST (%), As (g/t), Sb (g/t)
- **Rows**: 4 composites (Green-Breccia, White-Breccia, Basalt, Argillized)
- **Cu Range**: 146 - 3961 g/t (0.0146% - 0.3961%)
- **Status**: ✓ Clean, usable
- **Notes**: High Cu in Argillized composite

### Dataset 2: Unnamed Deposit (Lines 11-60)
- **Project**: Unknown (multiple deposits: COP, CR, TAJ, NAP, Diorite, Andesite)
- **Format**: Deposit/Composite level
- **Columns**: Deposit, Composite, Ag (g/t), Au (g/t), S (%), Cu (%), Pb (%), Zn (%), Mn (%)
- **Rows**: 49 samples
- **Cu Range**: 0.01% - 1.14%
- **Status**: ✓ Clean, usable
- **Notes**: Multiple deposit types, good Cu variation

### Dataset 3: Geological Domain (Lines 62-101)
- **Project**: Unknown VMS/LFZ project
- **Format**: Domain-based variability samples
- **Columns**: Geological Domain, Met Sample Name, Sample Type, Cu (%), Fe (%), S (%), Zn (%), Au (g/t), Ag (g/t), As (%)
- **Rows**: 24 samples
- **Cu Range**: 0.01% - 3.51%
- **Status**: ✓ Clean, usable
- **Notes**: VMS and LFZ domains, includes waste sample

### Dataset 4: Simple P/S Samples (Lines 103-110)
- **Project**: Unknown
- **Format**: Simple sample designations (P1-P5, S1-S2)
- **Columns**: Sample, Cu (%), Au (ppm), Fe (%), Mo (%), Stot (%), S= (%)
- **Rows**: 7 samples
- **Cu Range**: 0.26% - 1.11%
- **Status**: ✓ Clean, usable
- **Notes**: Includes Mo data

### Dataset 5: POR/RHY Composites (Lines 112-121)
- **Project**: Unknown (Porphyry/Rhyolite?)
- **Format**: Composite assays
- **Columns**: Composite, Cu, Fe, Au, As, S, C, CuOX, CuCN (all % or g/t)
- **Rows**: 8 composites
- **Cu Range**: 0.24% - 1.03%
- **Status**: ✓ Clean, usable
- **Notes**: Includes Cu oxidation states (CuOX, CuCN)

### Dataset 6: Cactus Mine Sequential Cu (Lines 123-128)
- **Project**: ASCU Cactus Mine (from NI 43-101 PFS)
- **Format**: Sequential copper analysis
- **Columns**: Sample, Test No., Total Cu (%), Ca (%), Sequential Copper Analysis (CuAS, CuCN, Res. Cu), Acid Soluble Copper (%)
- **Rows**: 3 samples
- **Cu Range**: 0.066% - 0.553%
- **Status**: ✓ Specialized format, but usable
- **Notes**: Focus on Cu speciation

### Dataset 7: DPM FSU Flotation Data (Lines 130-137)
- **Project**: DPM FSU (from technical report 2025-10-29)
- **Format**: Flotation product distribution
- **Columns**: Product, Weight (%), Assays (Cu, As, S, Au, Ag %), Distribution (%)
- **Rows**: 4 products (concentrate, scav tail, ro tail, combined tail)
- **Cu Range**: 0.025% - 1.68%
- **Status**: ⚠ Different format (flotation products not feed composition)
- **Notes**: This is metallurgical test results, not ore composition

### Dataset 8: Box Athona ICP Scan (Lines 139-160)
- **Project**: Unknown (from NI 43-101 PEA - FOR - Oct 20 2025)
- **Format**: ICP elemental scan
- **Columns**: Element, Box (g/t), Athona (g/t)
- **Rows**: 19 elements
- **Cu**: 64.4 and 8.8 g/t
- **Status**: ✓ Clean, 2 samples
- **Notes**: Comprehensive elemental analysis

### Dataset 9: FQM Cayeli (Lines 162-166)
- **Project**: FQM Cayeli (from tech report Oct 2025)
- **Format**: Ore type composites
- **Columns**: Sample, Cu (%), Zn (%), Au (g/t), Fe (%), Pb (%), Co (ppm), Cr (ppm), Mn (%), Ni (ppm), Ca (%)
- **Rows**: 3 ore types
- **Cu Range**: 0.83% - 4.21%
- **Status**: ✓ Clean, usable
- **Notes**: VMS-style mineralization

### Dataset 10: KCA Leach Test Data (Lines 168-298+)
- **Project**: Unknown (KCA tests from Tech_Report.pdf)
- **Format**: Metallurgical test results with process parameters
- **Columns**: Very complex - includes Sample No., Test No., Description, Head/Tail Cu (mg/kg), Recovery %, Days of Leach, P80 Size, NaCN consumption, Lime, Cement
- **Rows**: Multiple test conditions (~20+)
- **Cu Range**: Various (in mg/kg format)
- **Status**: ⚠ SPECIAL - This contains BOTH composition AND process parameters!
- **Notes**: This is extremely valuable - includes process outputs with composition inputs

---

## STATUS SUMMARY (Composition Data)

**Examining first 300 lines reveals:**
- Clean compositional datasets: 6-7
- Specialized/process-mixed datasets: 2-3
- Unusable/incomplete datasets: 1

**Next steps**: Continue cataloging remaining ~1000 lines...

