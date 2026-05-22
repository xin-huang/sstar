from pathlib import Path

from sstar.configs.preprocess_config import PreprocessConfig


def test_preprocess_config_normalizes_output_file(tmp_path):
    rel_out = tmp_path / "nested" / "sample.features"
    cfg = PreprocessConfig(
        vcf_file=tmp_path / "input.vcf",
        chr_name="chr1",
        ref_ind_file=tmp_path / "ref.txt",
        tgt_ind_file=tmp_path / "tgt.txt",
        output_file=rel_out,
        win_len=1000,
        win_step=500,
        feature_config_file=tmp_path / "feature.yml",
    )

    assert cfg.output_file == Path(rel_out).expanduser().resolve()
