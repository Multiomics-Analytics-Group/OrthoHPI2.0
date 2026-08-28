"""
Exports a Zotero collection to a BibTeX file for the manuscript.

This is the local stand-in for Overleaf's Zotero integration: it pulls a
collection straight out of the desktop client's database, so the bibliography
is whatever Zotero says it is and no reference is ever typed twice.

The database is opened read-only through an `immutable=1` URI, which takes no
lock and writes no journal file, so this is safe to run while Zotero is open.

    .venv/bin/python scripts/zotero_to_bib.py --collection postdoc \
        --out paper/refs.bib

Better BibTeX (retorque.re/zotero-better-bibtex) does the same job from inside
Zotero and rewrites the file on every edit rather than on demand. If it gets
installed, export the collection as Better BibLaTeX with "keep updated" ticked
and this script becomes redundant -- the citation keys it generates follow the
same author-year-word convention, so the .tex does not have to change.
"""
import argparse
import os
import re
import sqlite3
import sys
import unicodedata

DEFAULT_DB = os.path.expanduser('~/Zotero/zotero.sqlite')

# Zotero item type -> BibTeX entry type. Everything not named here falls back
# to @misc, which bibtex renders rather than drops.
#
# These are the classic BibTeX types, not biblatex's richer set: tectonic runs
# bibtex itself but has to shell out to biber, and its bundled biblatex (3.17)
# does not match the biber Homebrew ships (2.22). Classic bibtex also happens
# to be what most journal classes expect. So preprints and datasets become
# @misc rather than @online and @dataset.
ENTRY_TYPES = {
    'journalArticle': 'article',
    'preprint': 'misc',
    'dataset': 'misc',
    'book': 'book',
    'bookSection': 'incollection',
    'conferencePaper': 'inproceedings',
    'thesis': 'phdthesis',
    'report': 'techreport',
    'webpage': 'misc',
    'computerProgram': 'misc',
}

# Zotero field -> BibTeX field. Fields absent here are dropped: accessDate,
# libraryCatalog, rights and friends are catalogue metadata, not bibliography
# data, and only make the file noisier.
FIELDS = {
    'title': 'title',
    'publicationTitle': 'journal',
    'volume': 'volume',
    'issue': 'number',
    'pages': 'pages',
    'DOI': 'doi',
    'ISSN': 'issn',
    'url': 'url',
}

# Words too weak to identify a paper in a citation key.
STOPWORDS = {'a', 'an', 'the', 'of', 'on', 'in', 'for', 'and', 'to', 'with'}


def ascii_slug(text):
    """Strips accents and punctuation, the way a citation key needs."""
    text = unicodedata.normalize('NFKD', text)
    text = ''.join(c for c in text if not unicodedata.combining(c))

    return re.sub(r'[^A-Za-z0-9]', '', text).lower()


def read_collection(db_path, collection):
    """
    Returns the items of a named collection as
    (item_key, item_type, {field: value}, [surname, ...]).

    Both the field values and the creators come back in Zotero's own order,
    which for creators is author order -- `orderIndex` -- and matters.
    """
    uri = f'file:{db_path}?immutable=1'
    connection = sqlite3.connect(uri, uri=True)
    try:
        rows = connection.execute(
            'SELECT collectionID, collectionName FROM collections').fetchall()
        matches = [r for r in rows if r[1] == collection]
        if not matches:
            names = ', '.join(sorted(r[1] for r in rows))
            raise SystemExit(
                f'no collection named {collection!r} in {db_path}\n'
                f'available: {names}')
        collection_id = matches[0][0]

        items = connection.execute("""
            SELECT i.itemID, i.key, it.typeName
            FROM collectionItems ci
            JOIN items i ON ci.itemID = i.itemID
            JOIN itemTypes it ON i.itemTypeID = it.itemTypeID
            LEFT JOIN deletedItems di ON di.itemID = i.itemID
            WHERE ci.collectionID = ? AND di.itemID IS NULL
            ORDER BY ci.orderIndex, i.itemID
        """, (collection_id,)).fetchall()

        result = []
        for item_id, item_key, type_name in items:
            fields = dict(connection.execute("""
                SELECT f.fieldName, idv.value
                FROM itemData id
                JOIN fields f ON id.fieldID = f.fieldID
                JOIN itemDataValues idv ON id.valueID = idv.valueID
                WHERE id.itemID = ?
            """, (item_id,)).fetchall())

            creators = [row[0] for row in connection.execute("""
                SELECT CASE WHEN c.fieldMode = 1 THEN c.lastName
                            ELSE c.lastName || ', ' || c.firstName END
                FROM itemCreators ic
                JOIN creators c ON ic.creatorID = c.creatorID
                JOIN creatorTypes ct ON ic.creatorTypeID = ct.creatorTypeID
                WHERE ic.itemID = ? AND ct.creatorType = 'author'
                ORDER BY ic.orderIndex
            """, (item_id,)).fetchall()]

            result.append((item_key, type_name, fields, creators))

        return result
    finally:
        connection.close()


def year_of(fields):
    """Zotero stores dates as `2024-05-13 2024-05-13`, or just `2024`."""
    match = re.search(r'\b(1[6-9]\d{2}|20\d{2})\b', fields.get('date', ''))

    return match.group(1) if match else ''


def citation_key(fields, creators, taken):
    """
    author + year + first meaningful title word, e.g. `meitil2024analysis`.

    The same convention Better BibTeX's default produces, so swapping to it
    later does not invalidate the \\cite commands already written.
    """
    surname = ascii_slug(creators[0].split(',')[0]) if creators else 'anon'
    words = [w for w in re.findall(r'[A-Za-z]+', fields.get('title', ''))
             if w.lower() not in STOPWORDS]
    stem = f"{surname}{year_of(fields)}{ascii_slug(words[0]) if words else ''}"

    key, suffix = stem, ord('a')
    while key in taken:
        key = f'{stem}{chr(suffix)}'
        suffix += 1
    taken.add(key)

    return key


def escape(value):
    """
    Escapes the TeX specials that appear in real bibliographic data.

    Backslash first, or it would escape the backslashes it just introduced.
    Braces are left alone: Zotero users put them round protected capitals.
    """
    for char in '\\&%$#_':
        value = value.replace(char, '\\' + char)

    return ' '.join(value.split())


def format_entry(key, item_type, fields, creators):
    entry_type = ENTRY_TYPES.get(item_type, 'misc')

    pairs = []
    if creators:
        pairs.append(('author', ' and '.join(creators)))
    for zotero_field, bib_field in FIELDS.items():
        if fields.get(zotero_field):
            pairs.append((bib_field, fields[zotero_field]))
    if year_of(fields):
        pairs.append(('year', year_of(fields)))
    # A preprint has no journal, so say where it sits instead -- otherwise the
    # entry renders as a bare title and year with no way to find it.
    if item_type == 'preprint' and fields.get('repository'):
        pairs.append(('howpublished', f"Preprint, {fields['repository']}"))

    width = max(len(name) for name, _ in pairs)
    lines = [f'@{entry_type}{{{key},']
    lines += [f'  {name:<{width}} = {{{escape(value)}}},'
              for name, value in pairs]
    lines.append('}')

    return '\n'.join(lines)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--db', default=DEFAULT_DB,
                        help=f'Zotero database (default: {DEFAULT_DB})')
    parser.add_argument('--collection', default='postdoc',
                        help='collection to export (default: postdoc)')
    parser.add_argument('--out', default='paper/refs.bib',
                        help='output file (default: paper/refs.bib)')
    args = parser.parse_args()

    if not os.path.exists(args.db):
        raise SystemExit(f'no Zotero database at {args.db}')

    items = read_collection(args.db, args.collection)
    if not items:
        raise SystemExit(f'collection {args.collection!r} is empty')

    taken = set()
    entries = []
    for _, item_type, fields, creators in items:
        key = citation_key(fields, creators, taken)
        entries.append(format_entry(key, item_type, fields, creators))

    header = (f'% Generated by scripts/zotero_to_bib.py from the Zotero\n'
              f'% collection "{args.collection}". Do not edit: add or correct\n'
              f'% the reference in Zotero and run `make refs`.\n')

    os.makedirs(os.path.dirname(os.path.abspath(args.out)), exist_ok=True)
    with open(args.out, 'w') as handle:
        handle.write(header + '\n' + '\n\n'.join(entries) + '\n')

    print(f'{len(entries)} references -> {args.out}', file=sys.stderr)
    for key in sorted(taken):
        print(f'  {key}', file=sys.stderr)


if __name__ == '__main__':
    main()
