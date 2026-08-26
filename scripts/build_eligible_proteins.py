"""
Write eligible_proteins.parquet for a data directory that predates it.

The file holds the proteins the pipeline's filters passed, which is the background the
app tests a network's enrichment against. A full pipeline run writes it, but rebuilding a
data directory only for this file would repeat the orthology transfer, which takes hours
and would rewrite predictions.parquet from today's sources. This runs the filter stage
alone, off the same cached downloads, and writes nothing else.

    .venv/bin/python scripts/build_eligible_proteins.py [--config config.yml] [--data-dir data]

The filters read the secretome predictions and the DeepLoc output under the data
directory, so a directory whose predictions were built with other cut-offs than the ones
in pipeline/main.py will get a pool that does not match them.
"""
import argparse
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from pipeline import main


def check_inputs(data_dir):
    """
    Refuse a data directory that does not hold the filters' own inputs.

    A missing secretome file raises on its own, but a missing DeepLoc directory only
    warns and leaves the host proteins unfiltered, which would write a pool several times
    too large -- and a pool too large is read as a background, not as an error. The
    snapshot directories are the case that matters: they hold parquet outputs alone.
    """
    missing = [d for d in (os.path.join(data_dir, 'secretome'),
                           os.path.join(data_dir, main.DEEPLOC_ACCURATE_DIR))
               if not os.path.isdir(d)]
    if missing:
        raise SystemExit(f"{data_dir} does not hold the inputs the filters read: "
                         + ', '.join(missing) + "\nThe pool has to be built where the "
                         "secretome and DeepLoc predictions are.")


def build(config_file, data_dir):
    check_inputs(data_dir)
    print(f"Getting proteins for {config_file}...")
    proteins = main.get_proteins(config_file)
    print(f"  {sum(len(v) for v in proteins.values())} proteins before filtering")

    print("Applying secretome/tissue/DeepLoc filters...")
    proteins, _ = main.filter_proteins(config_file=config_file, data_dir=data_dir,
                                       proteins=proteins)
    print(f"  {len(proteins)} proteins after filtering")

    output_file = os.path.join(data_dir, 'eligible_proteins.parquet')
    main.save_eligible_proteins(proteins=proteins, output_file=output_file)
    print(f"Wrote {output_file}")


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--config', default='config.yml', help='configuration file to read')
    parser.add_argument('--data-dir', default='data', help='directory to write into')
    args = parser.parse_args()

    build(config_file=args.config, data_dir=args.data_dir)
