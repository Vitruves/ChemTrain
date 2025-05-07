import os
import subprocess
import tempfile
import shutil
import csv
import pytest

@pytest.fixture(scope="module")
def chemtrain_bin():
    bin_path = os.path.abspath(os.path.join(os.path.dirname(__file__), "../bin/chemtrain"))
    assert os.path.isfile(bin_path) and os.access(bin_path, os.X_OK), f"chemtrain binary not found at {bin_path}"
    return bin_path

@pytest.fixture
def temp_csv_files():
    temp_dir = tempfile.mkdtemp()
    input_path = os.path.join(temp_dir, "input.csv")
    output_path = os.path.join(temp_dir, "output.csv")
    yield input_path, output_path
    shutil.rmtree(temp_dir)

def write_csv(path, rows):
    with open(path, "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerows(rows)

def read_csv(path):
    with open(path, "r", newline="") as f:
        return list(csv.reader(f))

def test_csv_roundtrip_simple(chemtrain_bin, temp_csv_files):
    input_path, output_path = temp_csv_files
    rows = [
        ["smiles", "name"],
        ["CCO", "ethanol"],
        ["CC(=O)O", "acetic acid"],
        ["C1=CC=CC=C1", "benzene"]
    ]
    write_csv(input_path, rows)
    result = subprocess.run([chemtrain_bin, "-i", input_path, "-o", output_path], capture_output=True, text=True)
    assert result.returncode == 0, f"chemtrain failed: {result.stderr}"
    assert os.path.isfile(output_path), "Output CSV not created"
    out_rows = read_csv(output_path)
    assert len(out_rows) >= 2, "Output CSV too short"
    assert out_rows[0][0].lower().startswith("smiles"), "Output header missing or incorrect"
    assert any("ethanol" in row for row in out_rows), "Ethanol row missing in output"

def test_csv_handles_empty_file(chemtrain_bin, temp_csv_files):
    input_path, output_path = temp_csv_files
    write_csv(input_path, [])
    result = subprocess.run([chemtrain_bin, "-i", input_path, "-o", output_path], capture_output=True, text=True)
    assert result.returncode == 0, f"chemtrain failed on empty file: {result.stderr}"
    assert os.path.isfile(output_path), "Output CSV not created for empty input"
    out_rows = read_csv(output_path)
    assert out_rows == [] or all(len(row) == 0 for row in out_rows), "Output for empty input should be empty"

def test_csv_handles_large_file(chemtrain_bin, temp_csv_files):
    input_path, output_path = temp_csv_files
    rows = [["smiles", "name"]]
    for i in range(1000):
        rows.append([f"CC{i}", f"mol{i}"])
    write_csv(input_path, rows)
    result = subprocess.run([chemtrain_bin, "-i", input_path, "-o", output_path], capture_output=True, text=True)
    assert result.returncode == 0, f"chemtrain failed on large file: {result.stderr}"
    out_rows = read_csv(output_path)
    assert len(out_rows) >= 1000, "Output CSV missing rows for large input"
    assert any("mol999" in row for row in out_rows), "Last row missing in output"

def test_csv_handles_special_characters(chemtrain_bin, temp_csv_files):
    input_path, output_path = temp_csv_files
    rows = [
        ["smiles", "name"],
        ['C1=CC(=O)C=CC1=O', 'p-benzoquinone, "special" chars'],
        ['C(C(=O)O)N', "glycine, α-amino acid"],
        ['C1CC1', "cyclopropane\nnewline"]
    ]
    write_csv(input_path, rows)
    result = subprocess.run([chemtrain_bin, "-i", input_path, "-o", output_path], capture_output=True, text=True)
    assert result.returncode == 0, f"chemtrain failed on special chars: {result.stderr}"
    out_rows = read_csv(output_path)
    assert any('p-benzoquinone' in row[1] for row in out_rows), "Special char row missing"
    assert any('glycine' in row[1] for row in out_rows), "Unicode row missing"
    assert any('cyclopropane' in row[1] for row in out_rows), "Newline row missing"

def test_csv_handles_missing_values(chemtrain_bin, temp_csv_files):
    input_path, output_path = temp_csv_files
    rows = [
        ["smiles", "name"],
        ["CCO", ""],
        ["", "unknown"],
        ["CCN", None]
    ]
    write_csv(input_path, rows)
    result = subprocess.run([chemtrain_bin, "-i", input_path, "-o", output_path], capture_output=True, text=True)
    assert result.returncode == 0, f"chemtrain failed on missing values: {result.stderr}"
    out_rows = read_csv(output_path)
    assert any(row[0] == "CCO" for row in out_rows), "CCO row missing"
    assert any(row[1] == "unknown" for row in out_rows if len(row) > 1), "unknown row missing"

def test_csv_handles_header_case_insensitivity(chemtrain_bin, temp_csv_files):
    input_path, output_path = temp_csv_files
    rows = [
        ["SMILES", "Name"],
        ["CCO", "ethanol"],
        ["CCN", "ethylamine"]
    ]
    write_csv(input_path, rows)
    result = subprocess.run([chemtrain_bin, "-i", input_path, "-o", output_path], capture_output=True, text=True)
    assert result.returncode == 0, f"chemtrain failed: {result.stderr}"
    out_rows = read_csv(output_path)
    assert any("ethanol" in row for row in out_rows), "Case-insensitive header failed"

def test_csv_handles_extra_columns(chemtrain_bin, temp_csv_files):
    input_path, output_path = temp_csv_files
    rows = [
        ["smiles", "name", "extra1", "extra2"],
        ["CCO", "ethanol", "foo", "bar"],
        ["CCN", "ethylamine", "baz", "qux"]
    ]
    write_csv(input_path, rows)
    result = subprocess.run([chemtrain_bin, "-i", input_path, "-o", output_path], capture_output=True, text=True)
    assert result.returncode == 0, f"chemtrain failed: {result.stderr}"
    out_rows = read_csv(output_path)
    assert any("foo" in row for row in out_rows), "Extra columns not preserved"

def test_csv_handles_duplicate_smiles(chemtrain_bin, temp_csv_files):
    input_path, output_path = temp_csv_files
    rows = [
        ["smiles", "name"],
        ["CCO", "ethanol1"],
        ["CCO", "ethanol2"]
    ]
    write_csv(input_path, rows)
    result = subprocess.run([chemtrain_bin, "-i", input_path, "-o", output_path], capture_output=True, text=True)
    assert result.returncode == 0, f"chemtrain failed: {result.stderr}"
    out_rows = read_csv(output_path)
    ethanol_rows = [row for row in out_rows if "CCO" in row[0]]
    assert len(ethanol_rows) >= 2, "Duplicate SMILES rows not preserved"

def test_csv_handles_only_header(chemtrain_bin, temp_csv_files):
    input_path, output_path = temp_csv_files
    rows = [["smiles", "name"]]
    write_csv(input_path, rows)
    result = subprocess.run([chemtrain_bin, "-i", input_path, "-o", output_path], capture_output=True, text=True)
    assert result.returncode == 0, f"chemtrain failed: {result.stderr}"
    out_rows = read_csv(output_path)
    assert len(out_rows) == 1, "Output should only have header row"

def test_csv_handles_whitespace_in_fields(chemtrain_bin, temp_csv_files):
    input_path, output_path = temp_csv_files
    rows = [
        ["smiles", "name"],
        [" CCO ", "  ethanol  "],
        ["CCN", " ethylamine "]
    ]
    write_csv(input_path, rows)
    result = subprocess.run([chemtrain_bin, "-i", input_path, "-o", output_path], capture_output=True, text=True)
    assert result.returncode == 0, f"chemtrain failed: {result.stderr}"
    out_rows = read_csv(output_path)
    assert any("ethanol" in "".join(row) for row in out_rows), "Whitespace not handled"

def test_csv_handles_tab_delimited(chemtrain_bin, temp_csv_files):
    input_path, output_path = temp_csv_files
    # Write a tab-delimited file
    with open(input_path, "w", newline="") as f:
        f.write("smiles\tname\nCCO\tethanol\nCCN\tethylamine\n")
    result = subprocess.run([chemtrain_bin, "-i", input_path, "-o", output_path], capture_output=True, text=True)
    # Accept either success or a clear error about delimiter
    assert result.returncode == 0 or "delimiter" in result.stderr.lower(), "Tab-delimited file not handled or error not clear"

def test_csv_handles_windows_newlines(chemtrain_bin, temp_csv_files):
    input_path, output_path = temp_csv_files
    rows = [
        ["smiles", "name"],
        ["CCO", "ethanol"],
        ["CCN", "ethylamine"]
    ]
    with open(input_path, "w", newline="") as f:
        writer = csv.writer(f, lineterminator="\r\n")
        writer.writerows(rows)
    result = subprocess.run([chemtrain_bin, "-i", input_path, "-o", output_path], capture_output=True, text=True)
    assert result.returncode == 0, f"chemtrain failed: {result.stderr}"
    out_rows = read_csv(output_path)
    assert any("ethanol" in row for row in out_rows), "Windows newlines not handled"

def test_csv_handles_utf8_bom(chemtrain_bin, temp_csv_files):
    input_path, output_path = temp_csv_files
    rows = [
        ["smiles", "name"],
        ["CCO", "ethanol"]
    ]
    with open(input_path, "w", encoding="utf-8-sig", newline="") as f:
        writer = csv.writer(f)
        writer.writerows(rows)
    result = subprocess.run([chemtrain_bin, "-i", input_path, "-o", output_path], capture_output=True, text=True)
    assert result.returncode == 0, f"chemtrain failed: {result.stderr}"
    out_rows = read_csv(output_path)
    assert any("ethanol" in row for row in out_rows), "UTF-8 BOM not handled"

def test_csv_handles_trailing_newlines(chemtrain_bin, temp_csv_files):
    input_path, output_path = temp_csv_files
    rows = [
        ["smiles", "name"],
        ["CCO", "ethanol"],
        ["CCN", "ethylamine"]
    ]
    with open(input_path, "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerows(rows)
        f.write("\n\n")
    result = subprocess.run([chemtrain_bin, "-i", input_path, "-o", output_path], capture_output=True, text=True)
    assert result.returncode == 0, f"chemtrain failed: {result.stderr}"
    out_rows = read_csv(output_path)
    assert any("ethanol" in row for row in out_rows), "Trailing newlines not handled"

def test_csv_handles_large_number_of_columns(chemtrain_bin, temp_csv_files):
    input_path, output_path = temp_csv_files
    header = ["smiles"] + [f"col{i}" for i in range(100)]
    row = ["CCO"] + [str(i) for i in range(100)]
    rows = [header, row]
    write_csv(input_path, rows)
    result = subprocess.run([chemtrain_bin, "-i", input_path, "-o", output_path], capture_output=True, text=True)
    assert result.returncode == 0, f"chemtrain failed: {result.stderr}"
    out_rows = read_csv(output_path)
    assert any("CCO" in row for row in out_rows), "Large number of columns not handled" 