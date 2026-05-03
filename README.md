# varlociraptor-inspect

Visual inspection of Varlociraptor VCF records.

## Usage

Paste your Varlociraptor VCF record into the text box and click **Inspect**.

## Direct URL Link Functionality

You can link directly to a visualization by encoding VCF data as URL parameters.

**Parameters:**
- `PROB_<event>` — Event probability in PHRED scale (e.g. `PROB_SOMATIC=2.5`)
- `AFD_<sample>` — Allele frequency distribution (e.g. `AFD_tumor=0.0=0.01,0.5=10.5`)
- `OBS_<sample>` — Observation string (e.g. `OBS_tumor=31Rv.p.+**..`)

**Example:**

http://localhost:8501/?PROB_SOMATIC=2.5&PROB_GERMLINE=4.6&AFD_tumor=0.0=0.01,0.1=10.5&OBS_tumor=31Rv.p.+**..

This will immediately render the visualizations without needing to paste any VCF data manually.