"""Build a plain-text description of a variant record, for use as LLM chat context."""

from varlociraptor_inspect.plotting import AFDData, OBSData, ProbData
from collections.abc import Mapping


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
        f"Sample {sample_name}: allele frequency distribution over "
        f"{len(dist_entries)} candidate frequencies."
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

    def summarize(entries, label):
        if not entries:
            return f"  {label}: 0 reads."
        total = sum(e.count for e in entries)
        high_mapq = sum(e.count for e in entries if e.mapq == "High MAPQ")
        odds_breakdown: dict[str, int] = {}
        for e in entries:
            odds_breakdown[e.posterior_odds] = (
                odds_breakdown.get(e.posterior_odds, 0) + e.count
            )
        odds_str = ", ".join(
            f"{v} {k}" for k, v in sorted(odds_breakdown.items(), key=lambda kv: -kv[1])
        )
        return (
            f"  {label}: {total} reads total, {high_mapq} with high mapping quality. "
            f"Evidence strength breakdown: {odds_str}."
        )

    lines = [f"Sample {sample_name}: read-level observations."]
    lines.append(summarize(obs.ref_observations, "Reference allele"))
    lines.append(summarize(obs.alt_observations, "Alternative allele"))
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
