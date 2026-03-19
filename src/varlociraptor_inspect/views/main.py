import streamlit as st
import pysam
import tempfile
import os
import re
from varlociraptor_inspect import plotting


def main_view():
    st.set_page_config(
        page_title="Varlociraptor Inspect",
    )
    st.title("Varlociraptor Inspect")
    st.text("Visual inspection of Varlociraptor VCF records.")

    # Check for URL parameters
    query_params = st.query_params

    # Build VCF from URL parameters if present
    if query_params:
        record_text = build_vcf_from_params(query_params)
        if record_text:
            st.info("Loaded data from URL parameters")
    else:
        record_text = None

    # Load record from text input (override URL params if user types)
    user_input = st.text_area(
        "Paste your Varlociraptor VCF record here (including header lines starting with #)",
        value=record_text or "",
        height=200,
    )

    if user_input:
        try:
            # Check if header is included
            if not user_input.startswith("##fileformat"):
                # Auto-generate header
                record_text = generate_header(user_input)
            else:
                record_text = user_input

            # Parse and display
            display_plots(record_text)

        except Exception as e:
            st.error(f"Error parsing VCF record: {str(e)}")


def build_vcf_from_params(params):
    """Build a VCF record from URL query parameters"""
    # Extract PROB_ fields
    prob_fields = {}
    afd_fields = {}
    obs_fields = {}

    for key, value in params.items():
        if key.startswith("PROB_"):
            prob_fields[key] = value
        elif key.startswith("AFD_"):
            sample_name = key.replace("AFD_", "")
            afd_fields[sample_name] = value
        elif key.startswith("OBS_"):
            sample_name = key.replace("OBS_", "")
            obs_fields[sample_name] = value

    if not prob_fields:
        return None

    # Get all sample names
    sample_names = set(list(afd_fields.keys()) + list(obs_fields.keys()))
    if not sample_names:
        return None

    # Build INFO field
    info_parts = [f"{k}={v}" for k, v in prob_fields.items()]
    info_field = ";".join(info_parts)

    # Build sample columns
    format_field = "AF:AFD:DP:OBS"
    sample_columns = []

    for sample in sorted(sample_names):
        af = "0.5"  # Default
        afd = afd_fields.get(sample, "0.0=0.01")
        dp = "100"  # Default
        obs = obs_fields.get(sample, ".")

        sample_columns.append(f"{af}:{afd}:{dp}:{obs}")

    # Build complete data line
    data_line = (
        f"chr1\t1000\t.\tA\tT\t.\t.\t{info_field}\t{format_field}\t"
        + "\t".join(sample_columns)
    )

    return data_line


def generate_header(data_line):
    """Generate VCF header from data line"""
    lines = data_line.strip().split("\n")

    # Find column header and data
    column_header = None
    actual_data = None

    for line in lines:
        if line.startswith("#CHROM"):
            column_header = line
        elif not line.startswith("#") and line.strip():
            actual_data = line
            break

    if not actual_data:
        return data_line

    fields = actual_data.split("\t")
    chrom = fields[0]
    pos = int(fields[1])
    info_field = fields[7] if len(fields) > 7 else ""

    # Extract PROB_ fields
    prob_fields = re.findall(r"PROB_(\w+)=", info_field)

    # Build header
    header_lines = [
        "##fileformat=VCFv4.2",
        f"##contig=<ID={chrom},length={pos + 1000}>",
    ]

    # Add PROB_ INFO fields
    for prob_field in prob_fields:
        header_lines.append(f"##INFO=<ID=PROB_{prob_field},Number=.,Type=Float>")

    # Add FORMAT fields
    format_fields = [
        "##FORMAT=<ID=DP,Number=1,Type=Integer>",
        "##FORMAT=<ID=AF,Number=1,Type=Float>",
        "##FORMAT=<ID=AFD,Number=.,Type=String>",
        "##FORMAT=<ID=OBS,Number=1,Type=String>",
        "##FORMAT=<ID=HINTS,Number=.,Type=String>",
    ]
    header_lines.extend(format_fields)

    # Generate column header if missing
    if not column_header:
        num_samples = len(fields) - 9
        sample_names = [f"sample{i + 1}" for i in range(num_samples)]
        column_header = (
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t"
            + "\t".join(sample_names)
        )

    return "\n".join(header_lines) + "\n" + column_header + "\n" + actual_data


def display_plots(record_text):
    """Parse VCF and display plots"""
    tmp_fd, tmp_path = tempfile.mkstemp(suffix=".vcf", text=True)
    try:
        with os.fdopen(tmp_fd, "w") as tmp:
            tmp.write(record_text)

        with pysam.VariantFile(tmp_path) as vcf:
            record = next(vcf)
            sample_names = list(record.samples.keys())

            st.success(
                f"Successfully parsed VCF record at {record.chrom}:{record.pos} with {len(sample_names)} sample(s)"
            )

            # Display Event Probabilities
            st.header("Event Probabilities")
            chart1 = plotting.visualize_event_probabilities(record)
            st.altair_chart(chart1, use_container_width=True)

            # Display plots for each sample
            for idx, sample_name in enumerate(sample_names, 1):
                st.divider()
                st.header(f"Sample {idx}: {sample_name}")

                st.subheader("Allele Frequency Distribution")
                chart2 = plotting.visualize_allele_frequency_distribution(
                    record, sample_name
                )
                st.altair_chart(chart2, use_container_width=True)

                st.subheader("Observations")
                chart3 = plotting.visualize_observations(record, sample_name)
                st.altair_chart(chart3, use_container_width=True)
    finally:
        if os.path.exists(tmp_path):
            os.unlink(tmp_path)
