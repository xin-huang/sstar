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
from typing import Dict, Any, List
from scipy.spatial import distance_matrix
from gaishi.registries.stat_registry import STAT_REGISTRY
from gaishi.stats import GenericStatistic


@STAT_REGISTRY.register("sstar")
class Sstar(GenericStatistic):
    """
    S* (Sstar) statistic.

    Computes haplotype-based S* scores per individual from a target genotype
    matrix and variant positions. This implementation removes singletons,
    builds a physical-distance matrix and a genotype-distance matrix (L1),
    then performs dynamic programming to find the best-scoring chain.
    """

    @staticmethod
    def compute(
        **kwargs,
    ) -> Dict[str, Any]:
        """
        Computes S* scores for each target sample (column).

        Parameters
        ----------
        **kwargs
            tgt_gts : np.ndarray
                Target genotype matrix of shape (n_sites, n_samples).
            pos : np.ndarray
                Genomic positions for the n_sites variants.
            match_bonus : int, default 5000
                Bonus added when genotype distance is 0 and physical distance ≥ 10.
            max_mismatch : int, default 5
                Maximum allowed genotype distance to still consider a pair.
            mismatch_penalty : int, default -10000
                Penalty applied when 0 < genotype distance ≤ max_mismatch.

        Returns
        -------
        dict
            {'sstar': list[float]} — per-sample S* scores.

        Raises
        ------
        ValueError
            If any required key is missing.
        """
        ref_gts, tgt_gts, pos = Sstar.require(kwargs, "ref_gts", "tgt_gts", "pos")
        match_bonus = kwargs.get("match_bonus", 5000)
        max_mismatch = kwargs.get("max_mismatch", 5)
        mismatch_penalty = kwargs.get("mismatch_penalty", -10000)

        variants_not_in_ref = np.sum(ref_gts, axis=1) == 0
        sub_tgt_gts = tgt_gts[variants_not_in_ref]
        sub_pos = pos[variants_not_in_ref]

        scores = Sstar._calc_sstar(
            tgt_gt=sub_tgt_gts,
            pos=sub_pos,
            match_bonus=match_bonus,
            max_mismatch=max_mismatch,
            mismatch_penalty=mismatch_penalty,
        )
        return {"sstar": scores}

    @staticmethod
    def _calc_sstar(
        tgt_gt: np.ndarray,
        pos: np.ndarray,
        match_bonus: int = 5000,
        max_mismatch: int = 5,
        mismatch_penalty: int = -10000,
    ) -> List[float]:
        """
        Core S* computation (ArchIE variant).

        Parameters
        ----------
        tgt_gt : np.ndarray
            (n_sites, n_samples) target genotype matrix.
        pos : np.ndarray
            Positions of length n_sites (same row order as tgt_gt).
        match_bonus : int, optional
            Bonus for pairs with gd==0 and distance ≥ 10. Default: 5000.
        max_mismatch : int, optional
            Maximum gd allowed for a mismatch pair. Default: 5.
        mismatch_penalty : int, optional
            Penalty for 0 < gd ≤ max_mismatch; pairs with gd > max_mismatch are invalid. Default: -10000.

        Returns
        -------
        list[float]
            Per-sample S* scores.
        """

        def _create_matrices(gt: np.ndarray, positions: np.ndarray, idx: int):
            # Focal individual's non-zero loci
            hap = gt[:, idx]
            mask = hap != 0
            pos_sub = positions[mask]

            # Genotype sub-matrix across all individuals at those loci
            geno_matrix = gt[mask]  # (n_pos, n_samples)

            # Remove singletons across all individuals
            keep = geno_matrix.sum(axis=1) != 1
            geno_matrix = geno_matrix[keep]
            pos_sub = pos_sub[keep]

            # Genotype distance (L1) across loci
            gd_mat = distance_matrix(geno_matrix, geno_matrix, p=1)  # (n_pos, n_pos)

            # Physical distance matrix (absolute)
            p = pos_sub.astype(np.float64)
            pos_mat = np.tile(p, (p.size, 1))
            pd_mat = np.abs(pos_mat.T - pos_mat)

            return pd_mat, gd_mat, pos_sub

        def _calc_ind_sstar(
            pd_matrix: np.ndarray,
            gd_matrix: np.ndarray,
            positions_sub: np.ndarray,
            match_bonus: int,
            max_mismatch: int,
            mismatch_penalty: int,
        ) -> float:
            # Initialize scoring matrix based on pair rules
            invalid = pd_matrix < 10
            scores = pd_matrix.copy()
            scores[invalid] = -np.inf
            scores[(~invalid) & (gd_matrix == 0)] += match_bonus
            scores[(~invalid) & (gd_matrix > 0) & (gd_matrix <= max_mismatch)] = (
                mismatch_penalty
            )
            scores[(~invalid) & (gd_matrix > max_mismatch)] = -np.inf

            n = positions_sub.size
            if n == 0:
                return 0.0

            max_scores = np.full(n, -np.inf, dtype=np.float64)
            # Each j, find best chain ending at j using pairs (i < j)
            for j in range(n):
                best = -np.inf
                for i in range(j):
                    best = max(best, max_scores[i] + scores[j, i], scores[j, i])
                max_scores[j] = best

            if np.isneginf(max_scores).all():
                return 0.0
            return float(np.nanmax(max_scores))

        _, n_samples = tgt_gt.shape
        sstar_scores: List[float] = []
        for i in range(n_samples):
            pdm, gdm, pos_sub = _create_matrices(tgt_gt, pos, i)
            score = _calc_ind_sstar(
                pdm, gdm, pos_sub, match_bonus, max_mismatch, mismatch_penalty
            )
            if score == -np.inf:
                score = 0.0
            sstar_scores.append(score)

        return sstar_scores
