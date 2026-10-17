"""Reference implementation of the helpers imported by 02_build_search_fasta.py.
Contract only — replace with your own module if you have one.
Exposes: read_fasta, write_fasta, dedup_by_sequence, tag
"""

def read_fasta(path):
    header, seq = None, []
    with open(path) as fh:
        for line in fh:
            line = line.rstrip("\n")
            if line.startswith(">"):
                if header is not None:
                    yield header, "".join(seq)
                header, seq = line[1:], []
            elif line:
                seq.append(line.strip())
    if header is not None:
        yield header, "".join(seq)

def write_fasta(records, path):
    with open(path, "w") as out:
        for header, seq in records:
            out.write(f">{header}\n{seq}\n")

def dedup_by_sequence(records):
    """Keep the FIRST occurrence of each sequence (published-before-proteome ordering
    preserves published provenance)."""
    seen = set()
    for header, seq in records:
        key = seq.upper()
        if key in seen:
            continue
        seen.add(key)
        yield header, seq

def tag(seq_id, grp="CANON", src="uniprot"):
    return f"{seq_id}|grp={grp}|src={src}"
