"""Print the structure of an .rds file without R.

Streams the gzipped R serialization and reports the object tree — slots,
classes, dimensions, factor levels — while discarding the bulk numeric
payloads instead of holding them in memory.

Usage:
    python rds_peek.py FILE.rds [--depth N] [--max-items N] [--values N]
"""

import argparse
import gzip
import struct
import sys

# SEXP type codes from R's serialize.c
NILSXP, SYMSXP, LISTSXP, CLOSXP, ENVSXP, PROMSXP, LANGSXP = 0, 1, 2, 3, 4, 5, 6
SPECIALSXP, BUILTINSXP, CHARSXP, LGLSXP = 7, 8, 9, 10
INTSXP, REALSXP, CPLXSXP, STRSXP, DOTSXP = 13, 14, 15, 16, 17
VECSXP, EXPRSXP, BCODESXP, EXTPTRSXP, WEAKREFSXP, RAWSXP, S4SXP = 19, 20, 21, 22, 23, 24, 25

REFSXP, NILVALUE_SXP, GLOBALENV_SXP = 255, 254, 253
UNBOUNDVALUE_SXP, MISSINGARG_SXP, BASENAMESPACE_SXP = 252, 251, 250
NAMESPACESXP, PACKAGESXP, PERSISTSXP = 249, 248, 247
CLASSREFSXP, GENERICREFSXP, BCREPDEF, BCREPREF = 246, 245, 244, 243
EMPTYENV_SXP, BASEENV_SXP = 242, 241
ATTRLANGSXP, ATTRLISTSXP, ALTREP_SXP = 240, 239, 238

TYPENAME = {
    NILSXP: "NULL", SYMSXP: "symbol", LISTSXP: "pairlist", CLOSXP: "closure",
    ENVSXP: "environment", LANGSXP: "call", CHARSXP: "char", LGLSXP: "logical",
    INTSXP: "integer", REALSXP: "double", CPLXSXP: "complex", STRSXP: "character",
    VECSXP: "list", EXPRSXP: "expression", RAWSXP: "raw", S4SXP: "S4",
}


class Node:
    """A parsed object: its R type, a short summary, and named children."""

    def __init__(self, rtype, summary="", children=None, tags=None):
        self.rtype = rtype
        self.summary = summary
        self.children = children if children is not None else []
        self.tags = tags if tags is not None else []
        self.attrs = {}


class CountingStream:
    """Wraps the decompressed stream so errors can report a byte offset."""

    def __init__(self, inner):
        self.inner = inner
        self.pos = 0

    def read(self, n):
        data = self.inner.read(n)
        self.pos += len(data)
        return data


class Reader:
    def __init__(self, stream, keep_values):
        self.f = CountingStream(stream)
        self.keep_values = keep_values
        self.refs = []
        self.trail = []

    # --- primitive reads -------------------------------------------------
    def int32(self):
        raw = self.f.read(4)
        if len(raw) < 4:
            raise EOFError(
                f"stream ended at byte {self.f.pos} "
                f"(trail: {' > '.join(self.trail[-12:])})"
            )
        return struct.unpack(">i", raw)[0]

    def length(self):
        n = self.int32()
        if n == -1:  # long vector: two more ints form a 64-bit length
            hi, lo = self.int32(), self.int32()
            return (hi << 32) | lo
        return n

    def skip(self, nbytes):
        remaining = nbytes
        while remaining > 0:
            chunk = self.f.read(min(remaining, 8 << 20))
            if not chunk:
                raise EOFError("truncated stream")
            remaining -= len(chunk)

    # --- header ----------------------------------------------------------
    def header(self):
        magic = self.f.read(2)
        if magic != b"X\n":
            raise ValueError(f"not an XDR-format RDS file (magic {magic!r})")
        version = self.int32()
        writer, minreader = self.int32(), self.int32()
        encoding = None
        if version >= 3:
            encoding = self.f.read(self.int32()).decode("ascii", "replace")

        def ver(v):
            return f"{v >> 16}.{(v >> 8) & 255}.{v & 255}"

        return version, ver(writer), ver(minreader), encoding

    # --- object dispatch -------------------------------------------------
    def read(self):
        flags = self.int32()
        rtype = flags & 0xFF
        has_attr = bool((flags >> 9) & 1)
        has_tag = bool((flags >> 10) & 1)
        levels = flags >> 12

        if rtype == NILVALUE_SXP or rtype == NILSXP:
            return Node(NILSXP, "NULL")
        if rtype == REFSXP:
            idx = flags >> 8
            if idx == 0:
                idx = self.int32()
            # Resolve to the referenced symbol's own name: these back-references
            # carry attribute and slot names, so a placeholder loses them.
            target = self.refs[idx - 1] if 0 < idx <= len(self.refs) else "?"
            return Node(SYMSXP, target)
        if rtype in (GLOBALENV_SXP, BASEENV_SXP, EMPTYENV_SXP, BASENAMESPACE_SXP):
            return Node(ENVSXP, "<env>")
        if rtype in (MISSINGARG_SXP, UNBOUNDVALUE_SXP):
            return Node(NILSXP, "<missing>")

        if rtype == SYMSXP:
            name = self.read().summary
            self.refs.append(name)
            return Node(SYMSXP, name)

        if rtype in (NAMESPACESXP, PACKAGESXP):
            self.refs.append("<namespace>")
            self.read()  # the name vector
            return Node(SYMSXP, "<namespace>")

        if rtype == CHARSXP:
            n = self.int32()
            if n == -1:
                return Node(CHARSXP, "NA")
            raw = self.f.read(n)
            return Node(CHARSXP, raw.decode("utf-8", "replace"))

        if rtype in (LISTSXP, LANGSXP, ATTRLISTSXP, ATTRLANGSXP, DOTSXP):
            return self.read_pairlist(has_attr, has_tag)

        if rtype == ENVSXP:
            self.refs.append("<env>")
            self.int32()  # locked
            for _ in range(3):  # enclos, frame, hashtab
                self.read()
            self.read()  # attributes
            return Node(ENVSXP, "<environment>")

        if rtype in (CLOSXP, PROMSXP):
            node = Node(rtype, "<function>" if rtype == CLOSXP else "<promise>")
            if has_attr:
                self.read()
            if has_tag:  # the closure environment, written only when non-NULL
                self.read()
            self.read()  # formals / promise value
            self.read()  # body / promise expression
            return node

        if rtype == ALTREP_SXP:
            info = self.read()
            state = self.read()
            self.read()  # attributes
            cls = info.children[0].summary if info.children else "?"
            # wrap_* ALTREPs hold the real vector as the car of their state
            # pairlist; unwrap it so callers see the data, not the wrapper.
            if state.rtype == LISTSXP and state.children:
                state = state.children[0]
            node = Node(state.rtype, f"{state.summary} [ALTREP {cls}]")
            node.attrs = state.attrs
            node.values = getattr(state, "values", None)
            node.length = getattr(state, "length", 0)
            return node

        if rtype == S4SXP:
            node = Node(S4SXP, "S4 object")
            if has_attr:
                node.attrs = self.read_attributes()
            return node

        if rtype == BCODESXP:
            self.read_bytecode()
            return Node(BCODESXP, "<bytecode>")

        if rtype in (EXTPTRSXP, WEAKREFSXP):
            self.refs.append("<extptr>")
            self.read()
            self.read()
            return Node(rtype, "<extptr>")

        if rtype == PERSISTSXP:
            self.read()
            return Node(NILSXP, "<persistent>")

        if rtype in (LGLSXP, INTSXP, REALSXP, CPLXSXP, RAWSXP, STRSXP, VECSXP, EXPRSXP):
            return self.read_vector(rtype, has_attr, levels)

        raise ValueError(f"unhandled SEXP type {rtype} (flags {flags:#x})")

    # --- byte code -------------------------------------------------------
    # Mirrors ReadBC/ReadBC1/ReadBCConsts/ReadBCLang in R's serialize.c. The
    # values are discarded, but the traversal must be exact or the stream
    # desynchronises for everything that follows.
    def read_bytecode(self):
        nreps = self.int32()
        self.bc_reps = [None] * (nreps + 1)
        return self.read_bc1()

    def read_bc1(self):
        self.read()  # the code, an INTSXP
        self.read_bc_consts()

    def read_bc_consts(self):
        n = self.int32()
        for _ in range(n):
            const_type = self.int32()
            if const_type == BCODESXP:
                self.read_bc1()
            elif const_type in (LANGSXP, LISTSXP, BCREPDEF, BCREPREF,
                                ATTRLANGSXP, ATTRLISTSXP):
                self.read_bc_lang(const_type)
            else:
                self.read()

    def read_bc_lang(self, const_type):
        if const_type == BCREPREF:
            self.int32()  # index into the reps table
            return
        if const_type in (BCREPDEF, LANGSXP, LISTSXP, ATTRLANGSXP, ATTRLISTSXP):
            has_attr = False
            if const_type == BCREPDEF:
                self.int32()  # position in the reps table
                const_type = self.int32()
            if const_type == ATTRLANGSXP:
                const_type, has_attr = LANGSXP, True
            elif const_type == ATTRLISTSXP:
                const_type, has_attr = LISTSXP, True
            if has_attr:
                self.read()  # attributes
            self.read()  # tag
            self.read_bc_lang(self.int32())  # car
            self.read_bc_lang(self.int32())  # cdr
            return
        self.read()

    def read_pairlist(self, has_attr, has_tag):
        node = Node(LISTSXP, "pairlist")
        if has_attr:
            node.attrs = self.read_attributes()
        tag = self.read().summary if has_tag else ""
        self.trail.append(tag or "<untagged>")
        car = self.read()
        self.trail.pop()
        node.children.append(car)
        node.tags.append(tag)
        cdr = self.read()
        if cdr.rtype == LISTSXP:
            node.children.extend(cdr.children)
            node.tags.extend(cdr.tags)
        return node

    def read_attributes(self):
        """The ATTRIB field is just another item: a tagged pairlist or NULL."""
        node = self.read()
        if node.rtype != LISTSXP:
            return {}
        return {
            tag: child
            for tag, child in zip(node.tags, node.children)
            if tag
        }

    def read_vector(self, rtype, has_attr, levels):
        n = self.length()
        node = Node(rtype)
        if rtype in (VECSXP, EXPRSXP):
            for i in range(n):
                self.trail.append(f"list[{i + 1}/{n}]")
                node.children.append(self.read())
                self.trail.pop()
            node.summary = f"list[{n}]"
        elif rtype == STRSXP:
            keep = min(n, self.keep_values)
            vals = [self.read().summary for _ in range(keep)]
            for _ in range(n - keep):
                self.read()
            node.summary = f"character[{n}]"
            node.values = vals
        else:
            width = {LGLSXP: 4, INTSXP: 4, REALSXP: 8, CPLXSXP: 16, RAWSXP: 1}[rtype]
            if n <= self.keep_values:
                raw = self.f.read(n * width)
                node.values = self.decode_numbers(rtype, raw, n)
            else:
                head = self.f.read(min(n, self.keep_values) * width)
                node.values = self.decode_numbers(rtype, head, min(n, self.keep_values))
                self.skip((n - min(n, self.keep_values)) * width)
            node.summary = f"{TYPENAME[rtype]}[{n}]"
        node.length = n
        if has_attr:
            node.attrs = self.read_attributes()
        return node

    @staticmethod
    def decode_numbers(rtype, raw, n):
        if rtype in (LGLSXP, INTSXP):
            return list(struct.unpack(f">{n}i", raw))
        if rtype == REALSXP:
            return list(struct.unpack(f">{n}d", raw))
        if rtype == RAWSXP:
            return list(raw)
        return ["<complex>"] * n


def describe(node, values_shown):
    """One-line description, enriched with class/dim attributes."""
    bits = []
    cls = node.attrs.get("class")
    if cls is not None and getattr(cls, "values", None):
        bits.append("/".join(str(v) for v in cls.values))
    elif node.rtype == S4SXP:
        bits.append("S4")
    else:
        bits.append(node.summary)

    dim = node.attrs.get("dim")
    if dim is not None and getattr(dim, "values", None):
        bits.append(" x ".join(str(v) for v in dim.values))
    elif node.rtype not in (S4SXP,) and cls is not None:
        bits.append(node.summary)

    vals = getattr(node, "values", None)
    if vals and values_shown:
        shown = ", ".join(str(v)[:30] for v in vals[:values_shown])
        more = " ..." if getattr(node, "length", 0) > len(vals[:values_shown]) else ""
        bits.append(f"[{shown}{more}]")

    lvl = node.attrs.get("levels")
    if lvl is not None and getattr(lvl, "values", None):
        shown = ", ".join(str(v) for v in lvl.values[:values_shown or 6])
        bits.append(f"levels: {shown} ...")

    return "  ".join(b for b in bits if b)


def walk(node, name, depth, out, max_depth, max_items, values_shown):
    indent = "  " * depth
    label = f"{name}: " if name else ""
    out.append(f"{indent}{label}{describe(node, values_shown)}")
    if depth >= max_depth:
        return

    named = list(node.attrs.items())
    # a plain list carries its element names in the `names` attribute
    names_attr = node.attrs.get("names")
    child_names = getattr(names_attr, "values", None) if names_attr else None
    for i, child in enumerate(node.children):
        tag = node.tags[i] if i < len(node.tags) and node.tags[i] else None
        if tag is None and child_names and i < len(child_names):
            tag = str(child_names[i])
        named.append((tag or f"[[{i + 1}]]", child))

    for i, (key, child) in enumerate(named):
        if key in ("names", "class", "dim"):
            continue
        if i >= max_items:
            out.append(f"{indent}  ... {len(named) - i} more")
            break
        walk(child, key, depth + 1, out, max_depth, max_items, values_shown)


def as_column(node):
    """Materialise a parsed vector as a plain Python list, decoding factors."""
    vals = getattr(node, "values", None)
    if vals is None:
        return None
    levels = node.attrs.get("levels")
    if levels is not None and getattr(levels, "values", None):
        lv = levels.values
        return [lv[v - 1] if 1 <= v <= len(lv) else None for v in vals]
    return list(vals)


def dump_tables(root, outdir):
    """Write the small, useful tables — cell annotations, genes, embeddings."""
    import csv
    import os

    os.makedirs(outdir, exist_ok=True)
    written = []

    def attr_of(node, name):
        return node.attrs.get(name)

    meta = attr_of(root, "meta.data")
    if meta is not None:
        names = getattr(attr_of(meta, "names"), "values", []) or []
        rows = getattr(attr_of(meta, "row.names"), "values", []) or []
        cols = [as_column(c) for c in meta.children]
        keep = [(n, c) for n, c in zip(names, cols) if c is not None]
        n_rows = len(rows) if rows else (len(keep[0][1]) if keep else 0)
        path = os.path.join(outdir, "meta_data.csv")
        with open(path, "w", newline="") as fh:
            w = csv.writer(fh)
            w.writerow(["cell"] + [n for n, _ in keep])
            for i in range(n_rows):
                w.writerow(
                    [rows[i] if i < len(rows) else i]
                    + [c[i] if i < len(c) else "" for _, c in keep]
                )
        written.append((path, f"{n_rows} cells x {len(keep)} columns"))

    ident = attr_of(root, "active.ident")
    if ident is not None:
        col = as_column(ident)
        cells = getattr(attr_of(ident, "names"), "values", []) or []
        path = os.path.join(outdir, "active_ident.csv")
        with open(path, "w", newline="") as fh:
            w = csv.writer(fh)
            w.writerow(["cell", "ident"])
            for i, v in enumerate(col or []):
                w.writerow([cells[i] if i < len(cells) else i, v])
        written.append((path, f"{len(col or [])} cells"))

    assays = attr_of(root, "assays")
    if assays is not None and assays.children:
        rna = assays.children[0]
        counts = attr_of(rna, "counts")
        if counts is not None:
            dn = attr_of(counts, "Dimnames")
            if dn is not None and dn.children:
                genes = getattr(dn.children[0], "values", []) or []
                path = os.path.join(outdir, "genes.txt")
                with open(path, "w") as fh:
                    fh.write("\n".join(str(g) for g in genes) + "\n")
                written.append((path, f"{len(genes)} genes"))

    reductions = attr_of(root, "reductions")
    if reductions is not None:
        rnames = getattr(attr_of(reductions, "names"), "values", []) or []
        for name, red in zip(rnames, reductions.children):
            emb = attr_of(red, "cell.embeddings")
            vals = getattr(emb, "values", None) if emb is not None else None
            if not vals:
                continue
            dim = attr_of(emb, "dim")
            dims = getattr(dim, "values", None) if dim is not None else None
            if not dims or len(dims) != 2:
                continue
            nr, nc = dims
            if nr * nc != len(vals):
                continue  # truncated during the scan; skip rather than mislead
            path = os.path.join(outdir, f"reduction_{name}.csv")
            with open(path, "w", newline="") as fh:
                w = csv.writer(fh)
                w.writerow([f"{name}_{j + 1}" for j in range(nc)])
                for i in range(nr):  # R stores matrices column-major
                    w.writerow([vals[j * nr + i] for j in range(nc)])
            written.append((path, f"{nr} x {nc}"))

    return written


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("path")
    ap.add_argument("--depth", type=int, default=4)
    ap.add_argument("--max-items", type=int, default=40)
    ap.add_argument("--values", type=int, default=6, help="values to print per vector")
    ap.add_argument("--dump", metavar="DIR", help="write annotation tables to DIR")
    ap.add_argument(
        "--keep",
        type=int,
        default=0,
        help="fully retain vectors up to this length (implied by --dump)",
    )
    args = ap.parse_args()

    opener = gzip.open
    with open(args.path, "rb") as probe:
        if probe.read(2) != b"\x1f\x8b":
            opener = open

    keep = args.keep or (1_000_000 if args.dump else 0)
    with opener(args.path, "rb") as stream:
        r = Reader(stream, keep_values=max(args.values, 12, keep))
        version, writer, minreader, encoding = r.header()
        print(f"RDS v{version}, written by R {writer}, needs R >= {minreader}, encoding {encoding}")
        print()
        try:
            root = r.read()
        except (EOFError, ValueError) as exc:
            print(f"PARSE ERROR: {exc}", file=sys.stderr)
            raise

    out = []
    walk(root, "", 0, out, args.depth, args.max_items, args.values)
    print("\n".join(out))

    if args.dump:
        print()
        for path, note in dump_tables(root, args.dump):
            print(f"wrote {path}  ({note})")


if __name__ == "__main__":
    sys.exit(main())
