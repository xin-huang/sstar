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


import numpy as np
from typing import Any


class Sstar:
    """
    S* (Sstar) statistic.

    Computes haplotype-based S* scores per individual from a target genotype
    matrix and variant positions. This implementation removes singletons,
    builds a physical-distance matrix and a genotype-distance matrix (L1),
    then performs dynamic programming to find the best-scoring chain.
    """

    @staticmethod
    def _require(kwargs: dict[str, Any], *names: str) -> tuple[Any, ...]:
        missing = [n for n in names if n not in kwargs]
        if missing:
            raise ValueError(f"Missing required arguments: {missing}")
        return tuple(kwargs[n] for n in names)

    @staticmethod
    def compute(
        **kwargs,
    ) -> dict[str, Any]:
        """
        Compute S* scores for each target sample.

        Parameters
        ----------
        **kwargs
            ref_gts : np.ndarray
                Reference genotype matrix of shape `(n_sites, n_ref_samples)`.
            tgt_gts : np.ndarray
                Target genotype matrix of shape `(n_sites, n_tgt_samples)`.
            pos : np.ndarray
                Genomic positions for the `n_sites` variants.
            match_bonus : int, optional
                Bonus added when genotype distance is 0 and physical distance is at
                least 10 bp. Default: 5000.
            max_mismatch : int, optional
                Maximum genotype distance allowed for a mismatch pair. Default: 5.
            mismatch_penalty : int, optional
                Penalty applied when genotype distance is greater than 0 and less
                than or equal to `max_mismatch`. Default: -10000.

        Returns
        -------
        dict
            Dictionary with:

            - `sstar`: per-sample S* scores.
            - `region_ind_SNP_number`: per-sample total SNP counts, calculated as
              the union of reference-polymorphic sites and focal target-carried
              reference-absent sites within the input region.

        Raises
        ------
        ValueError
            If any required key is missing.
        """
        ref_gts, tgt_gts, pos = Sstar._require(kwargs, "ref_gts", "tgt_gts", "pos")
        match_bonus = kwargs.get("match_bonus", 5000)
        max_mismatch = kwargs.get("max_mismatch", 5)
        mismatch_penalty = kwargs.get("mismatch_penalty", -10000)

        variants_in_ref = np.sum(ref_gts, axis=1) != 0
        variants_not_in_ref = np.sum(ref_gts, axis=1) == 0

        ref_pos = pos[variants_in_ref]
        tgt_pos = pos[variants_not_in_ref]
        tgt_gt = tgt_gts[variants_not_in_ref]

        scores, total_snps = Sstar._calc_sstar(
            ref_pos=ref_pos,
            tgt_pos=tgt_pos,
            tgt_gt=tgt_gt,
            match_bonus=match_bonus,
            max_mismatch=max_mismatch,
            mismatch_penalty=mismatch_penalty,
        )

        return {
            "sstar": scores,
            "region_ind_SNP_number": total_snps,
        }

    @staticmethod
    def _calc_sstar(
        ref_pos: np.ndarray,
        tgt_pos: np.ndarray,
        tgt_gt: np.ndarray,
        match_bonus: int = 5000,
        max_mismatch: int = 5,
        mismatch_penalty: int = -10000,
    ) -> tuple[list[float], list[int]]:
        """
        Core S* computation for one genomic region.

        Parameters
        ----------
        ref_pos : np.ndarray
            Positions of variants that are polymorphic in the reference population
            within the input region.
        tgt_pos : np.ndarray
            Positions of variants that are absent from the reference population
            within the input region. This array must have the same row order as
            `tgt_gt`.
        tgt_gt : np.ndarray
            Target genotype matrix of shape `(n_tgt_sites, n_tgt_samples)`,
            restricted to variants absent from the reference population.
        match_bonus : int, optional
            Bonus added to a pair when the focal target sample has genotype distance
            0 between the two variants and their physical distance is at least
            10 bp. Default: 5000.
        max_mismatch : int, optional
            Maximum genotype distance allowed for a mismatch pair. Default: 5.
        mismatch_penalty : int, optional
            Penalty applied when genotype distance is greater than 0 and less than
            or equal to `max_mismatch`; pairs with genotype distance greater than
            `max_mismatch` are invalid. Default: -10000.

        Returns
        -------
        tuple[list[float], list[int]]
            A tuple containing:

            - Per-sample S* scores.
            - Per-sample total SNP counts, where each count is
              `len(union(ref_pos, focal_tgt_pos))`.

            Here `focal_tgt_pos` is the subset of `tgt_pos` where the focal target
            sample has a non-zero genotype.
        """

        def _create_matrices(
            ref_positions: np.ndarray,
            tgt_positions: np.ndarray,
            gt: np.ndarray,
            idx: int,
        ) -> tuple[np.ndarray, np.ndarray, np.ndarray, int]:
            # Focal sample's carried reference-absent variants.
            hap = gt[:, idx]
            mask = hap != 0

            hap_sub = hap[mask]
            pos_sub = tgt_positions[mask]

            # Same meaning as the original `total_snps`:
            # union of ref-polymorphic sites and focal target-carried sites.
            total_snps = int(np.union1d(ref_positions, pos_sub).size)

            # Physical distance matrix.
            p = pos_sub.astype(np.float64)
            pd_mat = np.abs(p[:, None] - p[None, :])

            # Genotype distance for the focal sample only.
            h = hap_sub.astype(np.float64)
            gd_mat = np.abs(h[:, None] - h[None, :])

            return pd_mat, gd_mat, pos_sub, total_snps

        def _calc_ind_sstar(
            pd_matrix: np.ndarray,
            gd_matrix: np.ndarray,
            positions_sub: np.ndarray,
            match_bonus: int,
            max_mismatch: int,
            mismatch_penalty: int,
        ) -> float:
            n = positions_sub.size
            if n <= 2:
                return np.nan

            invalid = pd_matrix < 10

            scores = pd_matrix.copy()
            scores[invalid] = -np.inf
            scores[(~invalid) & (gd_matrix == 0)] += match_bonus
            scores[(~invalid) & (gd_matrix > 0) & (gd_matrix <= max_mismatch)] = (
                mismatch_penalty
            )
            scores[(~invalid) & (gd_matrix > max_mismatch)] = -np.inf

            max_scores = np.full(n, -np.inf, dtype=np.float64)

            # Each j: best S* chain ending at SNP j using pairs i < j.
            for j in range(n):
                best = -np.inf
                for i in range(j):
                    best = max(best, max_scores[i] + scores[j, i], scores[j, i])
                max_scores[j] = best

            if np.isneginf(max_scores).all():
                return np.nan

            return float(np.nanmax(max_scores))

        _, n_samples = tgt_gt.shape

        sstar_scores: list[float] = []
        total_snps_per_sample: list[int] = []

        for i in range(n_samples):
            pdm, gdm, pos_sub, total_snps = _create_matrices(
                ref_positions=ref_pos,
                tgt_positions=tgt_pos,
                gt=tgt_gt,
                idx=i,
            )

            score = _calc_ind_sstar(
                pd_matrix=pdm,
                gd_matrix=gdm,
                positions_sub=pos_sub,
                match_bonus=match_bonus,
                max_mismatch=max_mismatch,
                mismatch_penalty=mismatch_penalty,
            )

            sstar_scores.append(score)
            total_snps_per_sample.append(total_snps)

        return sstar_scores, total_snps_per_sample
