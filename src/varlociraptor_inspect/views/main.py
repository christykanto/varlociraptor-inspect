from itertools import chain
from typing import Sequence
import streamlit as st
import streamlit.components.v1 as components
import pysam
import tempfile
import os
import re
from varlociraptor_inspect import plotting
from varlociraptor_inspect.plotting import ProbData, AFDData, OBSData


def normalize_whitespace(text: str) -> str:
    """Normalize whitespace in VCF records - replace spaces with tabs in data lines."""
    lines = []
    for line in text.split("\n"):
        if line.startswith("#") or not line.strip():
            lines.append(line)
        else:
            lines.append(re.sub(r"[ \t]+", "\t", line.strip()))
    return "\n".join(lines)


def build_vcf_from_url_params() -> tuple[ProbData, list[AFDData], list[OBSData]] | None:
    """Build dataclass instances directly from URL query parameters."""
    params = st.query_params
    if not params:
        return None

    prob_fields: dict[str, str] = {}
    afd_fields: dict[str, str] = {}
    obs_fields: dict[str, str] = {}

    for key in params:
        if key.startswith("PROB_"):
            prob_fields[key] = params[key]
        elif key.startswith("AFD_"):
            afd_fields[key.removeprefix("AFD_")] = params[key]
        elif key.startswith("OBS_"):
            obs_fields[key.removeprefix("OBS_")] = params[key]

    if not prob_fields:
        return None

    prob_data = ProbData.from_dict(prob_fields)
    sample_names = sorted(set(chain(afd_fields.keys(), obs_fields.keys())))

    afd_data_list: list[AFDData] = []
    for sample in sample_names:
        afd = AFDData.from_string(sample, afd_fields.get(sample, ""))
        if afd is not None:
            afd_data_list.append(afd)

    obs_data_list: list[OBSData] = [
        OBSData.from_string(sample, obs_fields.get(sample, ""))
        for sample in sample_names
    ]

    return prob_data, afd_data_list, obs_data_list


def build_query_string(
    prob_fields: dict[str, str], afd_fields: dict[str, str], obs_fields: dict[str, str]
) -> str:
    """Build a URL query string from raw PROB/AFD/OBS field values."""
    from urllib.parse import urlencode

    params: list[tuple[str, str]] = []
    for event, phred in prob_fields.items():
        params.append((f"PROB_{event}", str(phred)))
    for sample, value in obs_fields.items():
        if value and value != ".":
            params.append((f"OBS_{sample}", str(value)))
    for sample, value in afd_fields.items():
        if value and value != ".":
            params.append((f"AFD_{sample}", str(value)))
    return urlencode(params)


def render_copy_link_button(query_string: str, key: str):
    """Render a button that copies a shareable link to this record to the clipboard."""
    components.html(
        f"""
        <button id="copy-link-btn-{key}" style="
            padding: 0.4rem 0.8rem;
            border-radius: 0.4rem;
            border: 1px solid rgba(49, 51, 63, 0.2);
            background-color: #f0f2f6;
            cursor: pointer;
            font-size: 0.9rem;
        ">📋 Copy link to this record</button>
        <span id="copy-link-status-{key}" style="margin-left: 0.5rem; font-size: 0.85rem;"></span>
        <script>
            const btn = document.getElementById("copy-link-btn-{key}");
            const status = document.getElementById("copy-link-status-{key}");
            btn.addEventListener("click", async () => {{
                const baseUrl = window.parent.location.origin + window.parent.location.pathname;
                const fullUrl = baseUrl + "?{query_string}";
                try {{
                    await navigator.clipboard.writeText(fullUrl);
                    status.textContent = "Copied!";
                    setTimeout(() => {{ status.textContent = ""; }}, 2000);
                }} catch (err) {{
                    status.textContent = "Copy failed";
                }}
            }});
        </script>
        """,
        height=45,
    )


def main_view():
    st.set_page_config(page_title="Varlociraptor Inspect")
    st.title("Varlociraptor Inspect")
    st.text("Visual inspection of Varlociraptor VCF records.")

    url_data = build_vcf_from_url_params()

    if url_data is not None:
        prob_data, afd_data_list, obs_data_list = url_data
        st.info("Loaded data from URL parameters")

        params = st.query_params
        prob_fields = {
            k.removeprefix("PROB_"): v
            for k, v in params.items()
            if k.startswith("PROB_")
        }
        afd_fields = {
            k.removeprefix("AFD_"): v for k, v in params.items() if k.startswith("AFD_")
        }
        obs_fields = {
            k.removeprefix("OBS_"): v for k, v in params.items() if k.startswith("OBS_")
        }
        render_copy_link_button(
            build_query_string(prob_fields, afd_fields, obs_fields), key="url"
        )

        st.header("Event Probabilities")
        event_chart = plotting.visualize_event_probabilities(prob_data)
        if event_chart is not None:
            st.altair_chart(event_chart, use_container_width=True)
        else:
            st.warning("No event probability data available.")

        if not afd_data_list and not obs_data_list:
            st.warning(
                "No sample data found in URL parameters. Showing event probabilities only."
            )
        else:
            afd_by_sample = {d.sample_name: d for d in afd_data_list}
            obs_by_sample = {d.sample_name: d for d in obs_data_list}
            all_sample_names = sorted(
                set(chain(afd_by_sample.keys(), obs_by_sample.keys()))
            )

            for idx, sample_name in enumerate(all_sample_names, 1):
                st.divider()
                st.header(f"Sample {idx}: {sample_name}")

                afd = afd_by_sample.get(sample_name)
                st.subheader("Allele Frequency Distribution")
                if afd is not None:
                    st.altair_chart(
                        plotting.visualize_allele_frequency_distribution(afd),
                        use_container_width=True,
                    )
                else:
                    st.warning(
                        f"No allele frequency data available for sample {sample_name}."
                    )

                obs = obs_by_sample.get(sample_name)
                st.subheader("Observations")
                if obs is not None:
                    st.altair_chart(
                        plotting.visualize_observations(obs), use_container_width=True
                    )
                else:
                    st.warning(
                        f"No observation data available for sample {sample_name}."
                    )

    else:
        record_text = st.text_area(
            "Paste your Varlociraptor VCF record here (including header lines starting with #)",
            value="",
            height=200,
        )

        if record_text:
            try:
                record_text = normalize_whitespace(record_text)

                if not record_text.startswith("##fileformat"):
                    lines = record_text.strip().split("\n")
                    column_header = None
                    data_line = None

                    for line in lines:
                        if line.startswith("#CHROM"):
                            column_header = line
                        elif not line.startswith("#") and line.strip():
                            data_line = line
                            break

                    if data_line:
                        fields = data_line.split("\t")

                        if len(fields) < 8:
                            raise ValueError(
                                "VCF record must have at least 8 tab-separated columns"
                            )

                        chrom = fields[0]
                        pos = int(fields[1])
                        info_field = fields[7] if len(fields) > 7 else ""
                        prob_fields = re.findall(r"PROB_(\w+)=", info_field)

                        header_lines = [
                            "##fileformat=VCFv4.2",
                            f"##contig=<ID={chrom},length={pos + 1000}>",
                        ]
                        for prob_field in prob_fields:
                            header_lines.append(
                                f"##INFO=<ID=PROB_{prob_field},Number=.,Type=Float>"
                            )
                        header_lines.extend(
                            [
                                "##FORMAT=<ID=DP,Number=1,Type=Integer>",
                                "##FORMAT=<ID=AF,Number=1,Type=Float>",
                                "##FORMAT=<ID=AFD,Number=.,Type=String>",
                                "##FORMAT=<ID=OBS,Number=1,Type=String>",
                                "##FORMAT=<ID=HINTS,Number=.,Type=String>",
                            ]
                        )

                        if not column_header:
                            num_samples = len(fields) - 9
                            sample_names = [
                                f"sample{i + 1}" for i in range(num_samples)
                            ]
                            column_header = (
                                "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t"
                                + "\t".join(sample_names)
                            )

                        record_text = (
                            "\n".join(header_lines)
                            + "\n"
                            + column_header
                            + "\n"
                            + data_line
                        )

                tmp_fd, tmp_path = tempfile.mkstemp(suffix=".vcf", text=True)
                try:
                    with os.fdopen(tmp_fd, "w") as tmp:
                        tmp.write(record_text)

                    with pysam.VariantFile(tmp_path) as vcf:
                        record = next(vcf)
                        sample_names = list(record.samples.keys())

                        st.success(
                            f"Successfully parsed VCF record at {record.chrom}:{record.pos} "
                            f"with {len(sample_names)} sample(s)"
                        )

                        prob_fields = {}
                        for info_key, info_value in record.info.items():
                            if not info_key.startswith("PROB_"):
                                continue
                            if isinstance(info_value, Sequence) and not isinstance(
                                info_value, str
                            ):
                                info_value = info_value[0] if info_value else None
                            if info_value is not None:
                                prob_fields[info_key.removeprefix("PROB_")] = str(
                                    info_value
                                )

                        afd_fields = {}
                        obs_fields = {}
                        for sname in sample_names:
                            s = record.samples[sname]
                            afd = s.get("AFD")
                            if isinstance(afd, Sequence) and not isinstance(afd, str):
                                afd = afd[0] if afd else None
                            if afd is not None:
                                afd_fields[str(sname)] = str(afd)

                            obs = s.get("OBS")
                            if isinstance(obs, Sequence) and not isinstance(obs, str):
                                obs = obs[0] if obs else None
                            if obs is not None:
                                obs_fields[str(sname)] = str(obs)

                        render_copy_link_button(
                            build_query_string(prob_fields, afd_fields, obs_fields),
                            key="paste",
                        )

                        prob_data = ProbData.from_record(record)
                        st.header("Event Probabilities")
                        event_chart = plotting.visualize_event_probabilities(prob_data)
                        if event_chart is not None:
                            st.altair_chart(event_chart, use_container_width=True)
                        else:
                            st.warning("No event probability data available.")

                        if not sample_names:
                            st.warning(
                                "No sample data found. Showing event probabilities only."
                            )
                        else:
                            for idx, sample_name in enumerate(sample_names, 1):
                                st.divider()
                                st.header(f"Sample {idx}: {sample_name}")

                                afd = AFDData.from_record(record, str(sample_name))
                                st.subheader("Allele Frequency Distribution")
                                if afd is not None:
                                    st.altair_chart(
                                        plotting.visualize_allele_frequency_distribution(
                                            afd
                                        ),
                                        use_container_width=True,
                                    )
                                else:
                                    st.warning(
                                        f"No allele frequency data available for sample {sample_name}."
                                    )

                                obs = OBSData.from_record(record, str(sample_name))
                                st.subheader("Observations")
                                st.altair_chart(
                                    plotting.visualize_observations(obs),
                                    use_container_width=True,
                                )
                finally:
                    if os.path.exists(tmp_path):
                        os.unlink(tmp_path)

            except Exception as e:
                st.error(f"Error parsing VCF record: {str(e)}")
