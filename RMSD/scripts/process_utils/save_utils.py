import os
import matplotlib.pyplot as plt
import pandas as pd


def plot_and_save_rmsd(path_to_rmsd_dir: str,
                       output_directory: str = ".",
                       output_name: str = "rmsd.png"
                       ) -> None:
    rmsd_csvs_first_panel = ["rmsd_protein.csv", "rmsd_dna.csv", "all.csv"]
    rmsd_csvs_second_panel = ["rmsd_dna_inner.csv", "rmsd_dna_outer.csv"]

    SMALL_SIZE = 16
    MEDIUM_SIZE = 18
    BIGGER_SIZE = 20

    plt.rc('font', size=SMALL_SIZE)  # controls default text sizes
    plt.rc('axes', titlesize=SMALL_SIZE)  # fontsize of the axes title
    plt.rc('axes', labelsize=BIGGER_SIZE)  # fontsize of the x and y labels
    plt.rc('xtick', labelsize=SMALL_SIZE)  # fontsize of the tick labels
    plt.rc('ytick', labelsize=SMALL_SIZE)  # fontsize of the tick labels
    plt.rc('legend', fontsize=MEDIUM_SIZE)  # legend fontsize

    colors = ["blue", "red", "black", "green", "magenta"]
    fig, axs = plt.subplots(1, 2, figsize=(41 / 2.54, 12 / 2.54), dpi=300)

    labels = [r"sec. struct. C$\rm\alpha$",
              r"DNA",
              r"sec. struct. C$\rm\alpha$ & DNA",
              r"DNA inner turn",
              r"DNA outer turn"]

    axis_limits = None
    for ind, rmsd_csv in enumerate([*rmsd_csvs_first_panel, *rmsd_csvs_second_panel]):
        df_rmsd = pd.read_csv(os.path.join(path_to_rmsd_dir, rmsd_csv))
        if rmsd_csv in rmsd_csvs_first_panel:
            index = 0
        else:
            index = 1
        axs[index].plot(df_rmsd["time_ns"], df_rmsd["rmsd"],
                    color=colors[ind], linewidth=1.0, label=labels[ind])
        if axis_limits is None:
            axis_limits = [0, df_rmsd["time_ns"].values[-1] + 10, 0, 10]

    for index in [0, 1]:
        axs[index].set(xlabel="time, ns", ylabel=r"RMSD, ${\rm\AA}$")
        axs[index].set_xlim(axis_limits[:2])
        axs[index].set_ylim(axis_limits[2:])

    ax_right_1 = axs[0].twinx()
    ax_right_1.set_yticks(axs[0].get_yticks())
    ax_right_1.set_ylim(axs[0].get_ylim())
    ax_right_1.set_yticklabels([])
    ax_right_2 = axs[1].twinx()
    ax_right_2.set_yticks(axs[1].get_yticks())
    ax_right_2.set_ylim(axs[1].get_ylim())
    ax_right_2.set_yticklabels([])
    leg_1 = axs[0].legend(loc="upper left", bbox_to_anchor=(0.305, 0.97))
    leg_2 = axs[1].legend(loc="upper left", bbox_to_anchor=(0.45, 0.97))

    for line in leg_1.get_lines():
        line.set_linewidth(4.0)

    for line in leg_2.get_lines():
        line.set_linewidth(4.0)

    os.makedirs(output_directory, exist_ok=True)
    plt.savefig(os.path.join(output_directory, output_name), bbox_inches="tight")
    plt.close()
