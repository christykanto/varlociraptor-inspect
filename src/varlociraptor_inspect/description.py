"""Build a plain-text description of a variant record, for use as LLM chat context."""

from collections.abc import Mapping

from varlociraptor_inspect.plotting import AFDData, OBSData, ProbData


def _describe_prob_data(prob_data: ProbData) -> str:
    if not prob_data.entries:
        return "No event probabilities are available for this record."
    lines = [
        "Event probabilities (posterior probability that each event is true, 0-1 scale):"
    ]
    for entry in sorted(prob_data.entries, key=lambda e: -e.probability):
        lines.append(f"- {entry.event}: {entry.probability:.4f}")
    return "\n".join(lines)


def _describe_afd(afd: AFDData | None, sample_name: str) -> str:
    if afd is None:
        return f"Sample {sample_name}: no allele frequency distribution data available."
    ml = next((e for e in afd.entries if e.entry_type == "ML Estimate"), None)
    dist_entries = [e for e in afd.entries if e.entry_type != "ML Estimate"]
    lines = [
        (
            f"Sample {sample_name}: allele frequency distribution over "
            f"{len(dist_entries)} candidate frequencies."
        )
    ]
    if ml is not None:
        lines.append(
            f"  Maximum-likelihood estimate: allele frequency = {ml.allele_frequency:.3f}, "
            f"posterior probability = {ml.probability:.4f}."
        )
    return "\n".join(lines)


def _describe_obs(obs: OBSData | None, sample_name: str) -> str:
    if obs is None or (not obs.ref_observations and not obs.alt_observations):
        return f"Sample {sample_name}: no read observation data available."

    def category_breakdown(entries, field_name):
        counts: dict[str, int] = {}
        for e in entries:
            val = getattr(e, field_name)
            counts[val] = counts.get(val, 0) + e.count
        return ", ".join(
            f"{v} {k}" for k, v in sorted(counts.items(), key=lambda kv: -kv[1])
        )

    def summarize(entries, label):
        if not entries:
            return f"  {label}: 0 reads."
        total = sum(e.count for e in entries)
        odds_str = category_breakdown(entries, "posterior_odds")
        strand_str = category_breakdown(entries, "strand")
        return f"  {label}: {total} reads. Evidence strength: {odds_str}. Strand: {strand_str}."

    alt_total = sum(e.count for e in obs.alt_observations)
    ref_total = sum(e.count for e in obs.ref_observations)
    lines = [
        f"Sample {sample_name}: {alt_total} reads support the ALT (variant) allele, "
        f"{ref_total} reads support the REF/other allele "
        f"({alt_total + ref_total} total observations)."
    ]
    lines.append(
        f"IMPORTANT: the number of reads supporting THIS variant in sample "
        f"{sample_name} is exactly {alt_total}."
    )
    lines.append(summarize(obs.alt_observations, "ALT-supporting reads"))
    lines.append(summarize(obs.ref_observations, "REF-supporting reads"))
    return "\n".join(lines)


def build_record_description(
    prob_data: ProbData,
    afd_by_sample: Mapping[str, AFDData | None],
    obs_by_sample: Mapping[str, OBSData | None],
    variant_info: dict[str, str] | None = None,
) -> str:
    """Build a plain-text description of a variant record for use as LLM context."""
    sections = []
    if variant_info:
        sections.append(
            f"Variant: {variant_info.get('chrom', '?')}:{variant_info.get('pos', '?')} "
            f"{variant_info.get('ref', '?')}>{variant_info.get('alt', '?')}"
        )
    sections.append(_describe_prob_data(prob_data))

    sample_names = sorted(set(afd_by_sample) | set(obs_by_sample))
    for sample_name in sample_names:
        sections.append(_describe_afd(afd_by_sample.get(sample_name), sample_name))
        sections.append(_describe_obs(obs_by_sample.get(sample_name), sample_name))

    return "\n\n".join(sections)
