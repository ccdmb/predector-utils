"""
Stored data for use during runtime or testing.
Because dbCAN is an evolving database, we'll have to maintain several models
for each database release. The latest version will always we the default one.
"""

from pkg_resources import resource_filename

import json
from typing import List

import xgboost as xgb


def get_interesting_pfam_ids(version: str = "latest") -> List[str]:
    fname = resource_filename(__name__, "hmms_of_interest.json")
    with open(fname, "r") as handle:
        d = json.load(handle)

    assert isinstance(d, dict)

    assert "Pfam" in d

    if version == "latest":
        assert "latest" in d
        assert "Pfam" in d["latest"]
        version = d["latest"]["Pfam"]
        assert isinstance(version, str)

    pfids = d["Pfam"].get(version, None)
    if not isinstance(pfids, list):
        # TODO make this a proper error if we decide to actually support this.
        raise ValueError("Don't have a list of interesting ids for the "
                         "Pfam version specified.")

    assert all([isinstance(x, str) for x in pfids]), pfids

    return pfids


def get_interesting_dbcan_ids(version: str = "latest") -> List[str]:
    fname = resource_filename(__name__, "hmms_of_interest.json")
    with open(fname, "r") as handle:
        d = json.load(handle)

    assert isinstance(d, dict)

    assert "dbCAN" in d

    if version == "latest":
        assert "latest" in d
        assert "dbCAN" in d["latest"]
        version = d["latest"]["dbCAN"]
        assert isinstance(version, str)

    ids = d["dbCAN"].get(version, None)
    if not isinstance(ids, list):
        # TODO make this a proper error if we decide to actually support this.
        raise ValueError("Don't have a list of interesting ids for the "
                         "dbCAN version specified.")

    assert all([isinstance(x, str) for x in ids]), ids

    return ids


def get_ltr_model(load_config=False) -> xgb.Booster:
    model_fname = resource_filename(__name__, "learn_to_rank-model.json")
    model = xgb.Booster(model_file=model_fname)

    features_fname = resource_filename(__name__, "learn_to_rank-features.json")
    with open(features_fname, "r") as handle:
        features = json.load(handle)

    model.feature_names = features["feature_names"]
    model.feature_types = features["feature_types"]

    if load_config:
        config_fname = resource_filename(__name__, "learn_to_rank-config.json")
        with open(config_fname) as handle:
            config = handle.read()
        model.load_config(config)
    return model
