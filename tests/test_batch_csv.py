import os
import csv
import pytest
from polymerase_tm.batch import from_csv, to_csv

@pytest.fixture
def temp_csv(tmp_path):
    def _create_csv(filename, rows):
        path = tmp_path / filename
        with open(path, "w", newline="", encoding="utf-8") as f:
            if not rows:
                return str(path)
            fieldnames = list(rows[0].keys())
            writer = csv.DictWriter(f, fieldnames=fieldnames)
            writer.writeheader()
            for row in rows:
                writer.writerow(row)
        return str(path)
    return _create_csv

def test_from_csv_check_pair(temp_csv):
    path = temp_csv("test_pair.csv", [
        {"name": "p1", "fwd": "ATCGATCGATCG", "rev": "GCGCGCGCGCGC"},
        {"name": "p2", "fwd": "AATTCCGGAATT", "rev": "GTGCCTCCGAGCCAGCACC"}
    ])
    results = from_csv(path, action="check_pair")
    assert len(results) == 2
    assert results[0]["name"] == "p1"
    assert "ta" in results[0]
    assert results[1]["name"] == "p2"

def test_from_csv_tm(temp_csv):
    path = temp_csv("test_tm.csv", [
        {"name": "s1", "seq": "ATCGATCGATCG"},
        {"name": "s2", "seq": "GCGCGCGCGCGC"}
    ])
    results = from_csv(path, action="tm")
    assert len(results) == 2
    assert results[0]["tm"] == 47  # Q5 for ATCGATCGATCG
    assert results[1]["gc_pct"] == 100.0

def test_from_csv_protocol(temp_csv):
    path = temp_csv("test_prot.csv", [
        {"name": "p1", "fwd": "ATCGATCGATCG", "rev": "GCGCGCGCGCGC", "template": "ATCGATCGATCGGCGCGCGCGCGC"}
    ])
    results = from_csv(path, action="protocol")
    assert len(results) == 1
    assert "cycling" in results[0]
    assert results[0]["polymerase"] == "q5"

def test_from_csv_sdm(temp_csv):
    path = temp_csv("test_sdm.csv", [
        {"name": "mut1", "template": "ATGTCCCTGCTCTTCTCTCGATGCAA", "mutation": "M1A", "mode": "point"}
    ])
    results = from_csv(path, action="sdm")
    assert len(results) >= 1
    assert results[0]["name"] in ["mut1", "mut1_opt1"]
    assert results[0]["mutation"] == "M1A"
    assert "forward" in results[0]
    assert "reverse" in results[0]

def test_to_csv_flatten(temp_csv, tmp_path):
    data = [
        {"name": "test", "nested_list": [{"a": 1}, {"b": 2}], "nested_dict": {"x": 10}, "normal": "string"}
    ]
    out_path = str(tmp_path / "out.csv")
    to_csv(data, out_path)
    
    with open(out_path, "r", encoding="utf-8") as f:
        reader = list(csv.DictReader(f))
        assert len(reader) == 1
        assert reader[0]["normal"] == "string"
        assert "a=1; b=2" in reader[0]["nested_list"] or "{'a': 1}; {'b': 2}" in reader[0]["nested_list"]
        assert "x=10" in reader[0]["nested_dict"]
