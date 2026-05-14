from crypticip.config import load_config


def test_load_default():
    cfg = load_config()
    assert "paths" in cfg
    assert cfg.dotted("scoring.weights.depth") == 0.25


def test_load_with_organism():
    cfg = load_config(organism="yeast")
    assert cfg.dotted("organism") == "yeast"
    assert cfg.dotted("proteome.uniprot_id") == "UP000002311"


def test_load_with_validation_overlay():
    cfg = load_config("config/validation.yaml")
    assert "validation_set" in cfg
    assert cfg.dotted("validation_gate.adar2_ip6_pocket_rank_max") == 3


def test_config_hash_stable():
    a = load_config().hash()
    b = load_config().hash()
    assert a == b
    assert len(a) == 12
