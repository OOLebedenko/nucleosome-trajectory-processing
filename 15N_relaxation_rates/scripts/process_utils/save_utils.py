from pyxmolpp2 import calc_autocorr_order_2
from process_utils.fit import repeated_fit_auto_correlation, __multi_exp_f
from process_utils.plot import add_relpath_to_top_corner, settings_plot, get_autocorr_graph_label
from glob import glob
from tqdm import tqdm
import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages


def calc_and_save_acorr(path_to_vector_csv_files,
                        dt_ns,
                        acorr_func_limit=-1,
                        thumbling_time=None,
                        acorr_func=calc_autocorr_order_2,
                        out_dir="."):
    index = None
    for path_to_vector_file in tqdm(sorted(path_to_vector_csv_files)):

        out_name = os.path.basename(path_to_vector_file).split("_")[0] + ".csv"
        vectors_1 = pd.read_csv(path_to_vector_file).values
        acorr = acorr_func(vectors_1, limit=acorr_func_limit)

        if index is None:
            index = len(acorr)
            time_ns = np.linspace(0, index * dt_ns, index, endpoint=False)

        if (thumbling_time is not None) and (thumbling_time > 0):
            acorr *= np.exp(-time_ns / thumbling_time)

        os.makedirs(out_dir, exist_ok=True)
        pd.DataFrame({
            "time_ns": time_ns,
            "acorr": acorr
        }).to_csv(os.path.join(out_dir, out_name), index=False)


def fit_and_save_acorr_func(path_to_acorr_files,
                            bounds,
                            rname_list,
                            p0=None,
                            lag_spacing="log",
                            n_lag_points=None,
                            output_directory="./",
                            limit_ns=None):
    path_to_ccr_csv_files = sorted(glob(os.path.join(path_to_acorr_files, "*.csv")))
    for bound in bounds:
        tau_table = pd.DataFrame()
        for acorr_corr_file in tqdm(path_to_ccr_csv_files, desc=output_directory):
            df = pd.read_csv(acorr_corr_file)
            if limit_ns:
                df = df[df["time_ns"] <= limit_ns]

            time_ns, acorr = df["time_ns"].values, df["acorr"].values

            if lag_spacing == "log":
                lag_index = np.unique(
                    np.logspace(0, int(np.log10(time_ns.size)), n_lag_points, endpoint=False).astype(int))
                acorr = np.take(acorr, lag_index)
                time_ns = np.take(time_ns, lag_index)

            popt = repeated_fit_auto_correlation(acorr, time_ns, bound, p0)
            name = os.path.splitext(os.path.basename(acorr_corr_file))[0]
            amplitudes = popt[::2]
            taus = popt[1::2]
            order = (len(bound[0]) + 1) // 2

            rid = int(name.split("_")[0])

            popt_dict = {
                'rId': rid,
                'rName': rname_list[rid - 1],
                'limit_ns': limit_ns
            }

            popt_dict.update(
                {("exp-%d-a%d" % (order, i + 1)): a for i, a in enumerate(amplitudes)}
            )
            popt_dict.update(
                {("exp-%d-tau%d" % (order, i + 1)): tau for i, tau in enumerate(taus)}
            )

            tau_table = pd.concat([tau_table, pd.DataFrame(popt_dict, index=[0])])

        tau_table.sort_values(by=['rId'], inplace=True)
        os.makedirs(output_directory, exist_ok=True)
        tau_table.to_csv(os.path.join(output_directory, 'tau_%d_exp.csv' % order), index=False)


def plot_and_save_acorr_with_fit(path_to_fit_csv: str,
                                 path_to_acorr_csv: str,
                                 output_directory: str,
                                 ) -> None:
    os.makedirs(output_directory, exist_ok=True)
    exp_order_acorrs_csv = glob(os.path.join(path_to_fit_csv, "*.csv"))
    for tau_order_fit in sorted(exp_order_acorrs_csv):
        name = os.path.basename(tau_order_fit).split(".")[0]
        with PdfPages(os.path.join(output_directory, name + ".pdf")) as pdf:
            csv_fit = os.path.join(path_to_fit_csv, name + ".csv")
            fit = pd.read_csv(csv_fit)
            for (ind, (_, fit_line)) in enumerate(tqdm(fit.iterrows(), desc="plot"), 3):
                acorr_file = "{}/{:02d}.csv".format(path_to_acorr_csv, int(fit_line["rId"]))
                acorr_df = pd.read_csv(acorr_file)
                time, acorr = acorr_df["time_ns"], acorr_df["acorr"]
                graph_label = get_autocorr_graph_label(fit_line)
                fig, ax = settings_plot(graph_label)
                ax.set_title("Autocorrelation {rid} {rname}".format(
                    rid=fit_line["rId"],
                    rname=fit_line["rName"],
                ))

                # ax.set_xlim(-3, 400)
                ax.set_ylim(-0.1, 1.1)
                ax.grid(color="grey", alpha=0.3)

                amplitude = fit_line.filter(like='-a')
                tau = fit_line.filter(like='-tau')
                ax.plot(time, acorr)
                ax.plot(time, __multi_exp_f(time, amplitude,
                                            tau, C=0))

                ax.axvline(fit_line["limit_ns"], color="palegreen", ls="--")

                add_relpath_to_top_corner(fig)
                pdf.savefig(fig)
                plt.close(fig)
