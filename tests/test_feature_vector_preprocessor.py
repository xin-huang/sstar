import numpy as np

from sstar.feature_vector_preprocessor import FeatureVectorPreprocessor


def test_run_calls_sstar_compute_and_formats_output(tmp_path, monkeypatch):
    ref_ind = tmp_path / "ref.ind"
    tgt_ind = tmp_path / "tgt.ind"
    feat_yaml = tmp_path / "feature.yaml"

    ref_ind.write_text("ref1\n")
    tgt_ind.write_text("tgt1\n")
    feat_yaml.write_text(
        "sstar:\n"
        "  match_bonus: 123\n"
        "  max_mismatch: 4\n"
        "  mismatch_penalty: -7\n"
    )

    preprocessor = FeatureVectorPreprocessor(
        ref_ind_file=str(ref_ind),
        tgt_ind_file=str(tgt_ind),
        feature_config_file=str(feat_yaml),
    )

    captured = {}

    def fake_compute(**kwargs):
        captured.update(kwargs)
        return {"sstar": np.array([42.0])}

    monkeypatch.setattr("sstar.feature_vector_preprocessor.Sstar.compute", fake_compute)

    out = preprocessor.run(
        chr_name="chr1",
        start=10,
        end=20,
        ploidy=2,
        is_phased=False,
        ref_gts=np.array([[0], [1]]),
        tgt_gts=np.array([[1], [1]]),
        pos=np.array([11, 19]),
    )

    assert set(captured.keys()) == {
        "ref_gts",
        "tgt_gts",
        "pos",
        "match_bonus",
        "max_mismatch",
        "mismatch_penalty",
    }
    assert out == [
        {
            "Chromosome": "chr1",
            "Start": 10,
            "End": 20,
            "Sample": "tgt1",
            "Sstar": 42.0,
        }
    ]
