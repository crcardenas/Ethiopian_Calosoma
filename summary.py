import pandas as pd
import argparse
import seaborn as sns
import matplotlib.pyplot as plt
import plotly.express as px
import os

sns.set_theme()


def load_data(files: list) -> pd.DataFrame:
    # Load and concatenate multiple TSV files.
    dfs = []
    for file in files:
        df = pd.read_csv(file, sep='\t')
        df.columns = df.columns.str.strip().str.replace(' ', '_')
        df['source_file'] = os.path.basename(file)
        dfs.append(df)
    return pd.concat(dfs, ignore_index=True)


def identify_char_cols(df: pd.DataFrame) -> list:
    # Identify ambiguous character columns.
    return [c for c in ['n', '-', '?', 'N'] if c in df.columns]

def summarize_by(df: pd.DataFrame, group_col: str, char_cols: list) -> pd.DataFrame:
    # Summarize by taxa or loci.
    df = df.copy()
    required = ['Sequence_length', 'Missing_percent', 'AT_content', 'GC_content', 'Undetermined_characters']
    missing = [col for col in required if col not in df.columns]
    if missing:
        raise ValueError(f"Missing required columns in input: {missing}")

    df['Undetermined_mean'] = df[char_cols].sum(axis=1)
    summary = df.groupby(group_col).agg({
        'Sequence_length': 'mean',
        'Undetermined_characters': 'sum',
        'Undetermined_mean': 'mean',
        'Missing_percent': 'mean',
        'AT_content': 'mean',
        'GC_content': 'mean'
    })
    summary['count'] = df.groupby(group_col).size()
    summary.rename(columns={
        'Sequence_length': 'SeqLen_mean',
        'Missing_percent': 'Missing_percent_mean',
        'AT_content': 'AT_content_mean',
        'GC_content': 'GC_content_mean'
    }, inplace=True)
    return summary


def flag_high_missing(summary: pd.DataFrame, threshold: float) -> pd.DataFrame:
    return summary[summary['Missing_percent_mean'] >= threshold]


def plot_scatter_seaborn(x, y, xlabel, ylabel, title, threshold, filename, sizes=None):
    g = sns.relplot(x=x, y=y, size=sizes, kind='scatter', height=6, aspect=1.33)
    ax = g.ax if hasattr(g, 'ax') else g.axes[0][0]
    ax.axhline(y=threshold, linestyle='--', color='red', label=f'Threshold ({threshold}%)')
    ax.set(xlabel=xlabel, ylabel=ylabel, title=title)
    ax.legend()
    g.figure.savefig(filename, dpi=300, bbox_inches='tight')
    plt.close(g.figure)
    print(f"Saved plot: {filename}")


def plot_3d(summary: pd.DataFrame, output_html: str, group_col: str, color: str = None):
    df = summary.reset_index()
    fig = px.scatter_3d(
        df,
        x='SeqLen_mean',
        y='GC_content_mean',
        z='Missing_percent_mean',
        color=color if color in df.columns else None,
        hover_name=group_col,
        opacity=0.85,
        labels={
            'SeqLen_mean': 'Mean Seq Length',
            'GC_content_mean': 'Mean GC Content',
            'Missing_percent_mean': 'Mean Missing %'
        },
        title=f'{group_col} 3D: Length vs GC vs Missing %'
    )
    fig.write_html(output_html)
    print(f"Saved 3D plot: {output_html}")


def main():
    parser = argparse.ArgumentParser(
        description='Summarize and visualize alignment stats.'
        'requires "AMAS.py summary --by-taxon ..." to have been run on a directory of aligned loci before running this script.' )
    parser.add_argument('-i', '--input', required=True, nargs='+', help='Input TSV file(s)')
    parser.add_argument('-o', '--output', required=True, help='Output file prefix')
    parser.add_argument('--threshold', type=float, default=50, help='Missing %% threshold for plots')
    parser.add_argument('--length_threshold', type=float, default=50, help='AT-vs-Length threshold line')
    parser.add_argument('--pdf', action='store_true', help='Save plots as PDF instead of PNG')
    parser.add_argument('--sample_data', help='CSV file mapping taxon_name to sample type')
    parser.add_argument('--loci_3d', action='store_true', help='Generate loci 3D plot')
    parser.add_argument('--taxa_3d', action='store_true', help='Generate taxa 3D plot')
    parser.add_argument('--plot_2d', action='store_true', help='Generate 2D scatter plots')
    args = parser.parse_args()

    ext = 'pdf' if args.pdf else 'png'

    print("Loading summary files")
    df = load_data(args.input)
    print("Loaded columns:", df.columns.tolist())

    char_cols = identify_char_cols(df)
    loci = summarize_by(df, 'Alignment_name', char_cols)
    taxa = summarize_by(df, 'Taxon_name', char_cols)

    loci.to_csv(f"{args.output}_locus_summary.csv")
    taxa.to_csv(f"{args.output}_taxon_summary.csv")
    print("Saved summary CSVs")

    high_loci = flag_high_missing(loci, args.threshold)
    high_taxa = flag_high_missing(taxa, args.threshold)
    if not high_taxa.empty:
        print("High-missing taxa:", high_taxa.index.tolist())
    if not high_loci.empty:
        print("High-missing loci:", high_loci.index.tolist())

    if args.plot_2d:
        # Loci plots
        plot_scatter_seaborn(loci['SeqLen_mean'], loci['Missing_percent_mean'],
                             'Mean Alignment Length', 'Mean Missing %',
                             'Loci: Length vs Missing', args.threshold,
                             f"{args.output}_loci_len_missing.{ext}", sizes=loci['count'])

        plot_scatter_seaborn(loci['GC_content_mean'], loci['Missing_percent_mean'],
                             'Mean GC Content', 'Mean Missing %',
                             'Loci: GC vs Missing', args.threshold,
                             f"{args.output}_loci_gc_missing.{ext}", sizes=loci['count'])

        plot_scatter_seaborn(loci['AT_content_mean'], loci['Missing_percent_mean'],
                             'Mean AT Content', 'Mean Missing %',
                             'Loci: AT vs Missing', args.threshold,
                             f"{args.output}_loci_at_missing.{ext}", sizes=loci['count'])

        plot_scatter_seaborn(loci['GC_content_mean'], loci['SeqLen_mean'],
                             'Mean GC Content', 'Mean Alignment Length',
                             'Loci: GC vs Length', args.threshold,
                             f"{args.output}_loci_gc_len.{ext}", sizes=loci['count'])

        plot_scatter_seaborn(loci['AT_content_mean'], loci['SeqLen_mean'],
                             'Mean AT Content', 'Mean Alignment Length',
                             'Loci: AT vs Length', args.length_threshold,
                             f"{args.output}_loci_at_len.{ext}", sizes=loci['count'])

        # Taxa plot
        plot_scatter_seaborn(taxa['SeqLen_mean'], taxa['Missing_percent_mean'],
                             'Mean Seq Length', 'Mean Missing %',
                             'Taxa: Length vs Missing', args.threshold,
                             f"{args.output}_taxa_len_missing.{ext}", sizes=taxa['count'])

    # 3D plots
    if args.loci_3d:
        plot_3d(loci, f"{args.output}_loci_3d.html", 'alignment_name')

    if args.taxa_3d:
        color = None
        if args.sample_data:
            sample_map = pd.read_csv(args.sample_data, names=['Taxon_name', 'Type'])
            taxa = taxa.reset_index().merge(sample_map, on='Taxon_name', how='left').set_index('taxon_name')
            color = 'Type'
        plot_3d(taxa, f"{args.output}_taxa_3d.html", 'Taxon_name', color)


if __name__ == '__main__':
    main()