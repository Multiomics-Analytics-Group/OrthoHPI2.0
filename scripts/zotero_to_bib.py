'''Exports a Zotero collection to a BibTeX file for the manuscript.'''
import argparse
import os
import re
import sqlite3
import sys
import unicodedata

DEFAULT_DB = os.path.expanduser('~/Zotero/zotero.sqlite')

# Zotero item type -> BibTeX entry type; anything else becomes @misc. Classic BibTeX types,
# since tectonic's biblatex does not match Homebrew's biber
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

# Zotero field -> BibTeX field; fields absent here are dropped
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

# words too weak to identify a paper in a citation key
STOPWORDS = {'a', 'an', 'the', 'of', 'on', 'in', 'for', 'and', 'to', 'with'}


def ascii_slug(text):
    '''Strips accents and punctuation, the way a citation key needs.'''
    text = unicodedata.normalize('NFKD', text)
    text = ''.join(c for c in text if not unicodedata.combining(c))

    return re.sub(r'[^A-Za-z0-9]', '', text).lower()


def read_collection(db_path, collection):
    '''
    Returns the items of a named collection as (item_key, item_type, {field: value},
    [surname, ...]).
    '''
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
    '''Zotero stores dates as `2024-05-13 2024-05-13`, or just `2024`.'''
    match = re.search(r'\b(1[6-9]\d{2}|20\d{2})\b', fields.get('date', ''))

    return match.group(1) if match else ''


def citation_key(fields, creators, taken):
    '''author + year + first meaningful title word, e.g. `meitil2024analysis`.'''
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
    '''Escapes the TeX specials that appear in real bibliographic data.'''
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
    # a preprint has no journal, so say where it sits instead
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
