## Usage

Paste your Varlociraptor VCF record into the text box. The visualization renders automatically once a valid record is detected.

## Direct URL Link Functionality

You can link directly to a visualization by encoding VCF data as URL parameters.

**Parameters:**
- `PROB_<event>` — Event probability in PHRED scale (e.g. `PROB_SOMATIC=2.5`)
- `AFD_<sample>` — Allele frequency distribution (e.g. `AFD_tumor=0.0=0.01,0.5=10.5`)
- `OBS_<sample>` — Observation string (e.g. `OBS_tumor=31Rv.p.+**..`)

**Example:**

https://varlociraptor.github.io/varlociraptor-inspect/?PROB_SOMATIC=2.5&PROB_GERMLINE=4.6&AFD_tumor=0.0%3D0.01%2C0.1%3D10.5&OBS_tumor=31Rv.p.%2B%2A%2A..

This will immediately render the visualizations without needing to paste any VCF data manually. Note the parameter values above are URL-encoded (e.g. `+` becomes `%2B`) — the app's own "Copy link to this record" button generates these automatically, so you don't need to encode anything by hand.