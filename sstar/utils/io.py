# Copyright 2026 Xin Huang
#
# GNU General Public License v3.0
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program. If not, please see
#
#    https://www.gnu.org/licenses/gpl-3.0.en.html

import h5py
import multiprocessing
import numpy as np
import pandas as pd
from typing import Any, Mapping, Sequence, Union


def write_tsv(file_name: str, data_dict: dict, lock: multiprocessing.Lock) -> None:
    """
    Write the data dictionary to a TSV file.

    Parameters
    ----------
    file_name : str
        Path to the TSV file.
    data_dict : dict
        Dictionary containing the data to be written to the TSV file.
    lock : multiprocessing.Lock
        Lock for synchronizing multiprocessing operations.
    """
    converted_dict = {}
    for key, value in data_dict.items():
        if isinstance(value, np.ndarray):
            array_list = value.tolist()
            converted_dict[key] = array_list
        else:
            converted_dict[key] = value

    df = pd.DataFrame([converted_dict])

    with lock:
        with open(file_name, "a") as f:
            df.to_csv(f, header=False, index=False, sep="\t")


def write_h5(
    file_name: str,
    entries: Union[Mapping[str, Any], Sequence[Mapping[str, Any]]],
    ds_type: str,
    lock: multiprocessing.Lock,
) -> None:
    """
    Write inputs (and optional targets/coordinates) into an existing unified-schema HDF5 file.

    This function appends one or more entries into a pre-initialized HDF5 file opened in
    append/update mode by the caller (typically mode "a" inside the function). It writes
    only model inputs and, for `ds_type="train"`, the supervision targets. Prediction
    outputs such as logits are not written here.

    Data are stored in a unified, slice-friendly layout where the first dimension indexes
    replicates or windows. Sample order may differ between entries due to per-entry sorting.
    Therefore, per-entry row identity is stored as integer ids (`/index/ref_ids` and
    `/index/tgt_ids`) that point to global sample tables stored once in `/meta`.

    This function assumes the file has already been initialized (datasets created and
    `/meta` attributes present, including `n` and `n_written`). New entries are written
    starting at `row = /meta.attrs["n_written"]` and `n_written` is incremented accordingly.

    Parameters
    ----------
    file_name : str
        Output HDF5 path. The file must already exist and follow the unified schema.
    entries : Mapping[str, Any] or Sequence[Mapping[str, Any]]
        One entry or a list of entries. Each entry corresponds to one replicate/window.

        Required keys for both `ds_type="train"` and `ds_type="infer"`:
        - `Ref_genotype` : array_like, shape (N, L)
        - `Tgt_genotype` : array_like, shape (N, L)
        - `Gap_to_prev`  : array_like, shape (N, L)
        - `Gap_to_next`  : array_like, shape (N, L)
        - `Ref_sample`   : sequence of length N, sample names in row order
        - `Tgt_sample`   : sequence of length N, sample names in row order

        Additional required keys for `ds_type="train"`:
        - `Label`     : array_like, shape (N, L)
        - `Seed`      : scalar
        - `Replicate` : scalar

        Additional required keys for `ds_type="infer"`:
        - `Position`  : array_like, shape (L,)
    ds_type : {"train", "infer"}
        Dataset type. Controls whether training targets or inference coordinates are written.
    lock : multiprocessing.Lock
        Inter-process lock that serializes HDF5 writes.

    Notes
    -----
    Expected HDF5 schema (created elsewhere, e.g. by `initialize_h5`).

    Common datasets
    - `/data/Ref_genotype`  : uint32, shape (n, N, L)
    - `/data/Tgt_genotype`  : uint32, shape (n, N, L)
    - `/data/Gap_to_prev`   : int64,  shape (n, N, L)
    - `/data/Gap_to_next`   : int64,  shape (n, N, L)
    - `/index/ref_ids`      : uint32, shape (n, N)
    - `/index/tgt_ids`      : uint32, shape (n, N)
    - `/meta/ref_sample_table` : utf-8 strings, shape (K_ref,)
    - `/meta/tgt_sample_table` : utf-8 strings, shape (K_tgt,)
    - `/meta` attributes: `n`, `N`, `L`, `Chromosome`, `n_written`

    Training-only datasets (`ds_type="train"`)
    - `/targets/Label`      : uint8,  shape (n, N, L)
    - `/index/Seed`         : int64,  shape (n,)
    - `/index/Replicate`    : int64,  shape (n,)

    Inference-only datasets (`ds_type="infer"`)
    - `/coords/Position`    : int64,  shape (n, L)

    The integer row id arrays map each row in `Ref_genotype` and `Tgt_genotype` back to
    the global sample tables. This preserves per-entry sorting while still allowing efficient
    slicing across entries.
    """
    if ds_type not in ("train", "infer"):
        raise ValueError('ds_type must be "train" or "infer"')

    if isinstance(entries, Mapping):
        entries_list = [entries]
    else:
        entries_list = list(entries)
    if not entries_list:
        raise ValueError("entries is empty")

    e0 = entries_list[0]

    base_required = (
        "Ref_genotype",
        "Tgt_genotype",
        "Gap_to_prev",
        "Gap_to_next",
        "Ref_sample",
        "Tgt_sample",
    )
    for k in base_required:
        if k not in e0:
            raise KeyError(f"Missing key in entry[0]: {k}")

    if ds_type == "train":
        extra_required = ("Label", "Seed", "Replicate")
    else:
        extra_required = ("Position",)
    for k in extra_required:
        if k not in e0:
            raise KeyError(f"Missing key in entry[0]: {k}")

    ref0 = np.asarray(e0["Ref_genotype"])
    tgt0 = np.asarray(e0["Tgt_genotype"])
    if ref0.ndim != 2 or tgt0.ndim != 2:
        raise ValueError("Ref_genotype and Tgt_genotype must be 2D arrays (N, L)")
    if ref0.shape != tgt0.shape:
        raise ValueError(
            f"Ref_genotype shape {ref0.shape} != Tgt_genotype shape {tgt0.shape}"
        )
    N, L = ref0.shape
    n_new = len(entries_list)

    for i, e in enumerate(entries_list):
        ref = np.asarray(e["Ref_genotype"])
        tgt = np.asarray(e["Tgt_genotype"])
        if ref.shape != (N, L) or tgt.shape != (N, L):
            raise ValueError(f"Entry {i}: genotype shape mismatch, expected {(N, L)}")

        gp = np.asarray(e["Gap_to_prev"])
        gn = np.asarray(e["Gap_to_next"])
        if gp.shape != (N, L) or gn.shape != (N, L):
            raise ValueError(f"Entry {i}: gap shape mismatch, expected {(N, L)}")

        if len(e["Ref_sample"]) != N or len(e["Tgt_sample"]) != N:
            raise ValueError(
                f"Entry {i}: Ref_sample and Tgt_sample length must be N={N}"
            )

        if ds_type == "train":
            y = np.asarray(e["Label"])
            if y.shape != (N, L):
                raise ValueError(f"Entry {i}: Label shape mismatch, expected {(N, L)}")
        else:
            pos = np.asarray(e["Position"])
            if pos.ndim != 1 or pos.shape[0] != L:
                raise ValueError(f"Entry {i}: Position shape mismatch, expected {(L,)}")

    def _read_str_table(ds: h5py.Dataset) -> list[str]:
        arr = np.asarray(ds[()])
        out: list[str] = []
        for x in arr:
            if isinstance(x, (bytes, np.bytes_)):
                out.append(x.decode("utf-8"))
            else:
                out.append(str(x))
        return out

    with lock:
        with h5py.File(file_name, "a") as h5f:
            meta = h5f["/meta"]
            i0 = int(meta.attrs["n_written"])
            n_cap = int(meta.attrs["n"])
            if i0 + n_new > n_cap:
                raise ValueError(
                    f"HDF5 capacity exceeded: n_written={i0}, adding={n_new}, n={n_cap}"
                )

            ref_table = _read_str_table(h5f["/meta/ref_sample_table"])
            tgt_table = _read_str_table(h5f["/meta/tgt_sample_table"])
            ref_map = {name: idx for idx, name in enumerate(ref_table)}
            tgt_map = {name: idx for idx, name in enumerate(tgt_table)}

            for j, e in enumerate(entries_list):
                row = i0 + j

                h5f["/data/Ref_genotype"][row] = np.asarray(
                    e["Ref_genotype"], dtype=np.uint32
                )
                h5f["/data/Tgt_genotype"][row] = np.asarray(
                    e["Tgt_genotype"], dtype=np.uint32
                )
                h5f["/data/Gap_to_prev"][row] = np.asarray(
                    e["Gap_to_prev"], dtype=np.int64
                )
                h5f["/data/Gap_to_next"][row] = np.asarray(
                    e["Gap_to_next"], dtype=np.int64
                )

                h5f["/index/ref_ids"][row, :] = np.asarray(
                    [ref_map[str(s)] for s in e["Ref_sample"]], dtype=np.uint32
                )
                h5f["/index/tgt_ids"][row, :] = np.asarray(
                    [tgt_map[str(s)] for s in e["Tgt_sample"]], dtype=np.uint32
                )

                if ds_type == "train":
                    h5f["/targets/Label"][row] = np.asarray(e["Label"], dtype=np.uint8)
                    h5f["/index/Seed"][row] = int(e["Seed"])
                    h5f["/index/Replicate"][row] = int(e["Replicate"])
                else:
                    h5f["/coords/Position"][row] = np.asarray(
                        e["Position"], dtype=np.int64
                    )

            meta.attrs["n_written"] = i0 + n_new
            h5f.flush()


def initialize_h5(
    file_name: str,
    *,
    ds_type: str,
    num_genotype_matrices: int,
    N: int,
    L: int,
    chromosome: str,
    ref_table: list[str],
    tgt_table: list[str],
    compression: str = "lzf",
) -> None:
    """
    Initialize an HDF5 file in the unified schema.

    Parameters
    ----------
    file_name : str
        Path to the HDF5 file.
    ds_type : {"train", "infer"}
        Dataset type. "train" creates /targets/Label and /index/{Seed,Replicate}.
        "infer" creates /coords/Position.
    num_genotype_matrices : int
        Total number of genotype matrices (preallocated capacity).
    N : int
        Number of samples/rows per genotype matrix.
    L : int
        Number of sites/columns per genotype matrix.
    chromosome : str
        Chromosome identifier stored in /meta attrs.
    ref_table : list[str]
        Reference sample name table written once to /meta/ref_sample_table.
    tgt_table : list[str]
        Target sample name table written once to /meta/tgt_sample_table.
    compression : str, optional
        HDF5 dataset compression filter. Default: "lzf".

    Raises
    ------
    ValueError
        If ds_type is invalid or existing file metadata does not match inputs.
    """
    if ds_type not in ("train", "infer"):
        raise ValueError('ds_type must be "train" or "infer"')

    with h5py.File(file_name, "w") as h5f:
        h5f.require_group("/data")
        h5f.require_group("/index")
        meta = h5f.require_group("/meta")

        if ds_type == "train":
            h5f.require_group("/targets")
        else:
            h5f.require_group("/coords")

        # First-time initialization.
        meta.attrs["n"] = int(num_genotype_matrices)
        meta.attrs["N"] = int(N)
        meta.attrs["L"] = int(L)
        meta.attrs["Chromosome"] = str(chromosome)
        meta.attrs["n_written"] = 0

        str_dt = h5py.string_dtype(encoding="utf-8")
        h5f.create_dataset(
            "/meta/ref_sample_table",
            data=np.asarray(ref_table, dtype=object),
            dtype=str_dt,
        )
        h5f.create_dataset(
            "/meta/tgt_sample_table",
            data=np.asarray(tgt_table, dtype=object),
            dtype=str_dt,
        )

        h5f.create_dataset(
            "/data/Ref_genotype",
            shape=(num_genotype_matrices, N, L),
            dtype=np.uint32,
            chunks=(1, N, L),
            compression=compression,
        )
        h5f.create_dataset(
            "/data/Tgt_genotype",
            shape=(num_genotype_matrices, N, L),
            dtype=np.uint32,
            chunks=(1, N, L),
            compression=compression,
        )
        h5f.create_dataset(
            "/data/Gap_to_prev",
            shape=(num_genotype_matrices, N, L),
            dtype=np.int64,
            chunks=(1, N, L),
            compression=compression,
        )
        h5f.create_dataset(
            "/data/Gap_to_next",
            shape=(num_genotype_matrices, N, L),
            dtype=np.int64,
            chunks=(1, N, L),
            compression=compression,
        )

        h5f.create_dataset(
            "/index/ref_ids",
            shape=(num_genotype_matrices, N),
            dtype=np.uint32,
            chunks=(min(64, num_genotype_matrices), N),
            compression=compression,
        )
        h5f.create_dataset(
            "/index/tgt_ids",
            shape=(num_genotype_matrices, N),
            dtype=np.uint32,
            chunks=(min(64, num_genotype_matrices), N),
            compression=compression,
        )

        if ds_type == "train":
            h5f.create_dataset(
                "/targets/Label",
                shape=(num_genotype_matrices, N, L),
                dtype=np.uint8,
                chunks=(1, N, L),
                compression=compression,
            )
            h5f.create_dataset(
                "/index/Seed",
                shape=(num_genotype_matrices,),
                dtype=np.int64,
                chunks=(min(1024, num_genotype_matrices),),
                compression=compression,
            )
            h5f.create_dataset(
                "/index/Replicate",
                shape=(num_genotype_matrices,),
                dtype=np.int64,
                chunks=(min(1024, num_genotype_matrices),),
                compression=compression,
            )
        else:
            h5f.create_dataset(
                "/coords/Position",
                shape=(num_genotype_matrices, L),
                dtype=np.int64,
                chunks=(1, L),
                compression=compression,
            )
