from itertools import chain
from typing import Sequence
import streamlit as st
import streamlit.components.v1 as components
import pysam
import tempfile
import os
import re
import json
from varlociraptor_inspect import plotting
from varlociraptor_inspect.plotting import ProbData, AFDData, OBSData
from varlociraptor_inspect.description import build_record_description


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


def render_webllm_chat(description: str, key: str) -> None:
    """Render an in-browser chat widget (via WebLLM) that can answer questions about this record."""
    description_json = json.dumps(description)
    components.html(
        f"""
        <div id="webllm-chat-{key}" style="font-family: sans-serif; max-width: 700px;">
          <div id="webllm-status-{key}" style="font-size: 0.85rem; color: #555; margin-bottom: 6px;">
            Chat with this record using a small language model that runs entirely in your browser
            (nothing is sent to a server). Requires a WebGPU-capable browser (recent Chrome/Edge).
          </div>
          <div style="display:flex; gap: 8px; align-items:center; margin-bottom: 10px;">
            <select id="webllm-model-{key}" style="padding: 0.3rem; border-radius: 0.4rem;">
              <option value="Qwen2.5-0.5B-Instruct-q4f16_1-MLC">Qwen2.5 0.5B (fastest, ~0.4GB)</option>
              <option value="Llama-3.2-1B-Instruct-q4f16_1-MLC">Llama 3.2 1B (fast, ~0.9GB)</option>
              <option value="Phi-3.5-mini-instruct-q4f16_1-MLC">Phi-3.5 mini (better quality, ~2.2GB)</option>
            </select>
            <button id="webllm-load-btn-{key}" style="
                padding: 0.4rem 0.9rem;
                border-radius: 0.5rem;
                border: 1px solid rgba(49, 51, 63, 0.2);
                background-color: rgb(255, 255, 255);
                cursor: pointer;
            ">Load model &amp; start chat</button>
          </div>
          <div id="webllm-progress-{key}" style="font-size: 0.85rem; color: #888; margin-bottom: 10px;"></div>
          <div id="webllm-messages-{key}" style="
              display:none;
              max-height: 320px;
              overflow-y: auto;
              border: 1px solid rgba(49, 51, 63, 0.15);
              border-radius: 0.5rem;
              padding: 10px;
              margin-bottom: 8px;
              background: #fafafa;
          "></div>
          <div id="webllm-input-row-{key}" style="display:none; gap: 8px;">
            <input id="webllm-input-{key}" type="text" placeholder="Ask a question about this record..." style="
                flex: 1;
                padding: 0.5rem;
                border-radius: 0.4rem;
                border: 1px solid rgba(49, 51, 63, 0.3);
            "/>
            <button id="webllm-send-btn-{key}" style="
                padding: 0.4rem 0.9rem;
                border-radius: 0.5rem;
                border: 1px solid rgba(49, 51, 63, 0.2);
                background-color: rgb(255, 255, 255);
                cursor: pointer;
            ">Send</button>
          </div>
        </div>
        <script type="module">
          const statusEl = document.getElementById("webllm-status-{key}");
          const loadBtn = document.getElementById("webllm-load-btn-{key}");
          const modelSelect = document.getElementById("webllm-model-{key}");
          const progressEl = document.getElementById("webllm-progress-{key}");
          const messagesEl = document.getElementById("webllm-messages-{key}");
          const inputRow = document.getElementById("webllm-input-row-{key}");
          const inputEl = document.getElementById("webllm-input-{key}");
          const sendBtn = document.getElementById("webllm-send-btn-{key}");

          const systemPrompt = {description_json} + "\\n\\nYou are a helpful assistant answering questions about the variant record above. Always cite the specific numbers from the data above (probabilities, allele frequencies, read counts) in your answer. Do not give generic definitions - explain what these specific numbers imply about this specific variant. If something is not present in the data, say so instead of guessing.";
          let engine = null;
          let history = [{{ role: "system", content: systemPrompt }}];

          function appendMessage(role, text) {{
            const div = document.createElement("div");
            div.style.margin = "6px 0";
            div.style.whiteSpace = "pre-wrap";
            const label = document.createElement("strong");
            label.textContent = role === "user" ? "You: " : "Assistant: ";
            div.appendChild(label);
            const span = document.createElement("span");
            span.textContent = text;
            div.appendChild(span);
            messagesEl.appendChild(div);
            messagesEl.scrollTop = messagesEl.scrollHeight;
            return span;
          }}

          if (!navigator.gpu) {{
            statusEl.textContent = "Your browser does not support WebGPU, which is required for in-browser chat. Try a recent version of Chrome or Edge.";
            loadBtn.disabled = true;
          }}

          loadBtn.addEventListener("click", async () => {{
            loadBtn.disabled = true;
            modelSelect.disabled = true;
            try {{
              const webllm = await import("https://esm.run/@mlc-ai/web-llm");
              const modelId = modelSelect.value;
              progressEl.textContent = "Starting download...";
              engine = await webllm.CreateMLCEngine(modelId, {{
                initProgressCallback: (report) => {{
                  progressEl.textContent = report.text || "Loading model...";
                }},
              }});
              progressEl.textContent = "Model loaded. Ask a question below.";
              messagesEl.style.display = "block";
              inputRow.style.display = "flex";
            }} catch (err) {{
              progressEl.textContent = "Failed to load model: " + err.message;
              loadBtn.disabled = false;
              modelSelect.disabled = false;
            }}
          }});

          async function sendMessage() {{
            const text = inputEl.value.trim();
            if (!text || !engine) return;
            inputEl.value = "";
            sendBtn.disabled = true;
            appendMessage("user", text);
            history.push({{ role: "user", content: text }});
            const replySpan = appendMessage("assistant", "");
            try {{
              const stream = await engine.chat.completions.create({{
                messages: history,
                stream: true,
              }});
              let full = "";
              for await (const chunk of stream) {{
                const delta = chunk.choices[0]?.delta?.content || "";
                full += delta;
                replySpan.textContent = full;
                messagesEl.scrollTop = messagesEl.scrollHeight;
              }}
              history.push({{ role: "assistant", content: full }});
            }} catch (err) {{
              replySpan.textContent = "Error: " + err.message;
            }} finally {{
              sendBtn.disabled = false;
            }}
          }}

          sendBtn.addEventListener("click", sendMessage);
          inputEl.addEventListener("keydown", (e) => {{
            if (e.key === "Enter") {{
              sendMessage();
            }}
          }});
        </script>
        """,
        height=520,
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

            description = build_record_description(
                prob_data, afd_by_sample, obs_by_sample
            )
            st.divider()
            st.header("Chat with this record")
            render_webllm_chat(description, key="url")

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

                            afd_by_sample = {
                                str(s): AFDData.from_record(record, str(s))
                                for s in sample_names
                            }
                            obs_by_sample = {
                                str(s): OBSData.from_record(record, str(s))
                                for s in sample_names
                            }
                            description = build_record_description(
                                prob_data,
                                afd_by_sample,
                                obs_by_sample,
                                variant_info={
                                    "chrom": str(record.chrom),
                                    "pos": str(record.pos),
                                    "ref": str(record.ref),
                                    "alt": ",".join(str(a) for a in record.alts or ()),
                                },
                            )
                            st.divider()
                            st.header("Chat with this record")
                            render_webllm_chat(description, key="paste")
                finally:
                    if os.path.exists(tmp_path):
                        os.unlink(tmp_path)

            except Exception as e:
                st.error(f"Error parsing VCF record: {str(e)}")
