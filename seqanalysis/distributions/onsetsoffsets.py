import pathlib
import yaml
import re

import numpy as np
import polars as pl
import evfuncs
from IPython import embed
import matplotlib.pyplot as plt

from seqanalysis.util.logging import config_logging

log = config_logging()


def getAnnotations(files):
    """
    Get the annotations from the files

    Parameters
    ----------
    files : pathlib.Path
        Path to the files

    Returns
    -------
    polars_df : pl.DataFrame
        Dataframe with the annotations
    """
    labels = []
    onsets = []
    offsets = []
    filenames = []
    samplerates = []
    for f in files:
        notmat_dict = evfuncs.load_notmat(f)
        labels.extend(notmat_dict["labels"])
        onsets.extend(notmat_dict["onsets"])
        offsets.extend(notmat_dict["offsets"])
        filenames.extend([f.name] * len(notmat_dict["labels"]))
        samplerates.extend([notmat_dict["Fs"]])
    if np.unique(samplerates).shape[0] > 1:
        raise ValueError("The samplerates are not all the same")

    polars_df = pl.DataFrame(
        {
            "filename": filenames,
            "labels": labels,
            "onsets": onsets,
            "offsets": offsets,
        }
    )

    return polars_df


def generatePauses(polar_df, median: bool = True):
    """
    Generate pauses between syllables
    Parameters
    ----------
    polar_df : pl.DataFrame
        Dataframe with the labels
    median : bool, optional
        If True the median of the pauses is calculated for each syllable, by default True

    Returns
    -------
    pause_df : pl.DataFrame
        Dataframe with the pauses
    """
    files = polar_df["filename"].unique().sort().to_list()
    pause_df = pl.DataFrame(
        {
            "filename": [],
            "labels": [],
            "onsets": [],
            "offsets": [],
            "pauses_before": [],
            "pauses_after": [],
        },
        schema={
            "filename": str,
            "labels": str,
            "onsets": pl.Float64,
            "offsets": pl.Float64,
            "pauses_before": pl.Float64,
            "pauses_after": pl.Float64,
        },
    )

    for file in files:
        # get the rows of the dataframe where the filename is the current file
        file_df = polar_df.filter(pl.col("filename") == file)
        if file_df["labels"].unique().shape[0] == 1:
            log.error(f"File {file} has only one label skipping")
            continue
        # get all onsets and offsets from the dataframe
        onsets = file_df["onsets"].to_numpy()
        offsets = file_df["offsets"].to_numpy()

        # generate pauses after a syllable
        pauses_after = []
        for i in range(len(onsets) - 1):
            pauses_after.append(onsets[i + 1] - offsets[i])
        # add for the last syllable a nan value
        pauses_after.append(None)

        # generate pauses before
        pauses_before = []
        for i in range(1, len(onsets)):
            pauses_before.append(onsets[i] - offsets[i - 1])
        # add for the first syllable a nan value
        pauses_before.insert(0, None)

        # add the pauses to the dataframe
        file_df = file_df.with_columns(
            pl.Series("pauses_before", pauses_before),
            pl.Series("pauses_after", pauses_after),
        )

        # add the dataframe to the pause_df
        pause_df.extend(file_df)

    # adding median if needed!
    if median:
        # add mean values for each label
        median_values_df = pause_df.group_by(pl.col("labels")).agg(
            pl.col("pauses_before").median().alias("median_pause_before"),
            pl.col("pauses_after").median().alias("median_pause_after"),
        )
        # add the mean values to the sylor07ye05ause_d = pause_df.join(median_values_df, on="labels")

    return pause_df


def plot_pause(pause_df, syll_of_interest, mult_syllabels, pause_type="pauses_before"):
    pauses_onsets = {}
    fig, ax = plt.subplots()
    for pos, syllabels in enumerate(mult_syllabels):
        pauses_onsets[f"{syllabels}"] = []
        matches_idx = [
            (m.start(0), m.end(0))
            for m in re.finditer(syllabels, "".join(pause_df["labels"]))
        ]
        matches = re.findall(syllabels, "".join(pause_df["labels"]))
        matches_specific_syll = np.concatenate(
            [
                [m.start(0) for m in re.finditer(syll_of_interest, match)]
                for match in matches
            ]
        )
        log.info(f"Found {len(matches)} matches for {syllabels}")
        for match_specific_syl_idx, (indx_start, indx_stop) in enumerate(matches_idx):
            row_first_syll = pause_df.row(
                indx_start + matches_specific_syll[match_specific_syl_idx], named=True
            )
            if row_first_syll["labels"] != matches[0][matches_specific_syll[0]]:
                log.error("Label is not the same")
            pauses_onsets[f"{syllabels}"].append(row_first_syll[pause_type])

    ax.boxplot(
        pauses_onsets.values(),
        positions=np.arange(len(pauses_onsets)),
        patch_artist=True,
    )
    for i, (label, values) in enumerate(pauses_onsets.items()):
        x = np.random.normal(i, 0.01, size=len(values))
        ax.plot(x, values, "r.", alpha=0.5, label=f"{label} N:{len(values)}")
    ax.set_xticklabels(list(pauses_onsets.keys()))
    ax.set_title(f"{pause_type} k")
    ax.set_ylabel("Pause in ms")
    ax.legend()
    plt.show()


if __name__ == "__main__":
    with open(
        "/home/acfw/Projects/seqanalysis_py/seqanalysis/transition_diagram/gr88gr06.yaml"
    ) as f:
        cfg = yaml.load(f, Loader=yaml.FullLoader)
        f.close()
    folder_path = pathlib.Path(cfg["paths"]["folder_path"])
    files = list(folder_path.glob("**/[!syll*][!check]*.not.mat"))
    all_files_df = getAnnotations(files)
    onsets_offsets_df = generatePauses(all_files_df)
    plot_pause(onsets_offsets_df, "k", ["bnm+l+khtg", "bnm+l+kfg"], "pauses_before")
