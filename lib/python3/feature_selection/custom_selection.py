import warnings
import numpy as np
import rpy2.robjects as robjects
from rpy2.rinterface import RRuntimeWarning
from rpy2.robjects.packages import importr
from rpy2.robjects import numpy2ri
from sklearn.base import BaseEstimator
from sklearn.feature_selection.base import SelectorMixin
from sklearn.utils import check_array, check_X_y
from sklearn.utils.validation import check_is_fitted, check_memory
from .univariate_selection import BaseScorer

r_base = importr('base')
r_base.source('functions.R')
r_deseq2_vst_transform = robjects.globalenv['deseq2_vst_transform']
r_deseq2_feature_score = robjects.globalenv['deseq2_feature_score']
r_edger_filterbyexpr_mask = robjects.globalenv['edger_filterbyexpr_mask']
r_edger_tmm_logcpm_transform = robjects.globalenv['edger_tmm_logcpm_transform']
r_edger_feature_score = robjects.globalenv['edger_feature_score']
r_limma_voom_feature_score = robjects.globalenv['limma_voom_feature_score']
r_limma_feature_score = robjects.globalenv['limma_feature_score']
r_cfs_feature_idxs = robjects.globalenv['cfs_feature_idxs']
r_fcbf_feature_idxs = robjects.globalenv['fcbf_feature_idxs']
r_relieff_feature_score = robjects.globalenv['relieff_feature_score']
numpy2ri.activate()

def deseq2_vst_transform(X, y, geo_means, size_factors, disp_func, blind,
                         fit_type):
    if y is not None:
        xt, gm, sf, df = r_deseq2_vst_transform(
            X, y, blind=blind, fit_type=fit_type)
    else:
        xt, gm, sf, df = r_deseq2_vst_transform(
            X, geo_means=geo_means, size_factors=size_factors,
            disp_func=disp_func, blind=blind, fit_type=fit_type)
    return (np.array(xt, dtype=float), np.array(gm, dtype=float),
            np.array(sf, dtype=float), df)

def deseq2_feature_score(X, y, lfc, blind, fit_type):
    pv, pa, xt, gm, sf, df = r_deseq2_feature_score(X, y, lfc=lfc, blind=blind,
                                                    fit_type=fit_type)
    return (np.array(pv, dtype=float), np.array(pa, dtype=float),
            np.array(xt, dtype=float), np.array(gm, dtype=float),
            np.array(sf, dtype=float), df)

def edger_tmm_logcpm_transform(X, ref_sample, prior_count):
    xt, rs = r_edger_tmm_logcpm_transform(X, ref_sample=ref_sample,
                                          prior_count=prior_count)
    return np.array(xt, dtype=float), np.array(rs, dtype=float)

def edger_feature_score(X, y, lfc, robust, prior_count):
    pv, pa, xt, rs = r_edger_feature_score(X, y, lfc=lfc, robust=robust,
                                           prior_count=prior_count)
    return (np.array(pv, dtype=float), np.array(pa, dtype=float),
            np.array(xt, dtype=float), np.array(rs, dtype=float))

def limma_voom_feature_score(X, y, lfc, robust, prior_count):
    pv, pa, xt, rs = r_limma_voom_feature_score(X, y, lfc=lfc, robust=robust,
                                                prior_count=prior_count)
    return (np.array(pv, dtype=float), np.array(pa, dtype=float),
            np.array(xt, dtype=float), np.array(rs, dtype=float))

def fcbf_feature_idxs(X, y, threshold):
    idxs, scores = r_fcbf_feature_idxs(X, y, threshold=threshold)
    return np.array(idxs, dtype=int), np.array(scores, dtype=float)

def relieff_feature_score(X, y, num_neighbors, sample_size):
    return np.array(r_relieff_feature_score(X, y, num_neighbors=num_neighbors,
                                            sample_size=sample_size),
                    dtype=float)


class DESeq2(BaseEstimator, SelectorMixin):
    """DESeq2 feature selector and normalizer/transformer for RNA-seq data

    Parameters
    ----------
    k : int or "all" (default = "all")
        Number of top features to select. The "all" option bypasses selection,
        for use in a parameter search.

    lfc : float (default = 0)
        lfcThreshold absolute log fold-change minimum threshold.

    blind : bool (default = False)
        varianceStabilizingTransformation() blind option

    fit_type : str (default = local)
        estimateDispersions() fitType option

    memory : None, str or object with the joblib.Memory interface \
        (default = None)
        Used for internal caching. By default, no caching is done.
        If a string is given, it is the path to the caching directory.

    Attributes
    ----------
    pvals_ : array, shape (n_features,)
        Feature raw p-values.

    padjs_ : array, shape (n_features,)
        Feature adjusted p-values.

    geo_means_ : array, shape (n_features,)
        Feature geometric means.

    size_factors_ : array, shape (n_features,)
        RLE normalization size factors.

    disp_func_ : R/rpy2 function
        RLE normalization dispersion function.
    """
    def __init__(self, k='all', lfc=0, blind=False, fit_type='local',
                 memory=None):
        self.k = k
        self.lfc = lfc
        self.blind = blind
        self.fit_type = fit_type
        self.memory = memory

    def fit(self, X, y):
        """
        Parameters
        ----------
        X : array-like, shape = (n_samples, n_features)
            The training input counts data matrix.

        y : array-like, shape = (n_samples,)
            The target values (class labels in classification, real numbers in
            regression).

        Returns
        -------
        self : object
            Returns self.
        """
        X, y = check_X_y(X, y, dtype=np.int64)
        self._check_params(X, y)
        memory = check_memory(self.memory)
        (self.pvals_, self.padjs_, self._vst_data, self.geo_means_,
         self.size_factors_, self.disp_func_) = (
             memory.cache(deseq2_feature_score)(
                 X, y, lfc=self.lfc, blind=self.blind, fit_type=self.fit_type))
        return self

    def transform(self, X):
        """
        Parameters
        ----------
        X : array-like, shape = (n_samples, n_features)
            The input counts data matrix.

        Returns
        -------
        Xr : array of shape (n_samples, n_selected_features)
            DESeq2 median-of-ratios normalized VST transformed input data
            matrix with only the selected features.
        """
        check_is_fitted(self, '_vst_data')
        X = check_array(X, dtype=np.int64)
        if hasattr(self, '_train_done'):
            memory = check_memory(self.memory)
            X = memory.cache(deseq2_vst_transform)(
                X, y=None, geo_means=self.geo_means_,
                size_factors=self.size_factors_, disp_func=self.disp_func_,
                blind=self.blind, fit_type=self.fit_type)[0]
        else:
            X = self._vst_data
            self._train_done = True
        return super().transform(X)

    def inverse_transform(self, X):
        """
        Parameters
        ----------
        X : array-like, shape = (n_samples, n_features)
            The input data matrix.

        Returns
        -------
        Xr : array of shape (n_samples, n_original_features)
            `X` with columns of zeros inserted where features would have
            been removed by :meth:`transform`.
        """
        raise NotImplementedError("inverse_transform not implemented.")

    def _check_params(self, X, y):
        if not (self.k == "all" or 0 <= self.k <= X.shape[1]):
            raise ValueError(
                "k should be 0 <= k <= n_features; got %r."
                "Use k='all' to return all features."
                % self.k)

    def _get_support_mask(self):
        check_is_fitted(self, 'pvals_')
        mask = np.zeros_like(self.pvals_, dtype=bool)
        if self.k == 'all':
            mask = np.ones_like(self.pvals_, dtype=bool)
        elif self.k > 0:
            mask[np.argsort(self.pvals_, kind='mergesort')[:self.k]] = True
        return mask


class EdgeR(BaseEstimator, SelectorMixin):
    """edgeR feature selector and normalizer/transformer for RNA-seq data

    Parameters
    ----------
    k : int or "all" (default = "all")
        Number of top features to select. The "all" option bypasses selection,
        for use in a parameter search.

    lfc : float (default = 0)
        glmTreat absolute log fold-change minimum threshold.  Default value
        of 0 gives glmQLFTest results.

    robust : bool (default = True)
        estimateDisp and glmQLFit robust option

    prior_count : int (default = 1)
        Average count to add to each observation to avoid taking log of zero.
        Larger values for produce stronger moderation of the values for low
        counts and more shrinkage of the corresponding log fold changes.

    memory : None, str or object with the joblib.Memory interface \
        (default = None)
        Used for internal caching. By default, no caching is done.
        If a string is given, it is the path to the caching directory.

    Attributes
    ----------
    pvals_ : array, shape (n_features,)
        Feature raw p-values.

    padjs_ : array, shape (n_features,)
        Feature adjusted p-values.

    ref_sample_ : array, shape (n_features,)
        TMM normalization reference sample feature vector.
    """
    def __init__(self, k='all', lfc=0, robust=True, prior_count=1,
                 memory=None):
        self.k = k
        self.lfc = lfc
        self.robust = robust
        self.prior_count = prior_count
        self.memory = memory

    def fit(self, X, y):
        """
        Parameters
        ----------
        X : array-like, shape = (n_samples, n_features)
            The training input counts data matrix.

        y : array-like, shape = (n_samples,)
            The target values (class labels in classification, real numbers in
            regression).

        Returns
        -------
        self : object
            Returns self.
        """
        X, y = check_X_y(X, y, dtype=np.int64)
        self._check_params(X, y)
        memory = check_memory(self.memory)
        self.pvals_, self.padjs_, self._log_cpms, self.ref_sample_ = (
            memory.cache(edger_feature_score)(
                X, y, lfc=self.lfc, robust=self.robust,
                prior_count=self.prior_count))
        return self

    def transform(self, X):
        """
        Parameters
        ----------
        X : array-like, shape = (n_samples, n_features)
            The input counts data matrix.

        Returns
        -------
        X_r : array of shape (n_samples, n_selected_features)
            edgeR TMM normalized log-CPM transformed input data matrix with
            only the selected features.
        """
        check_is_fitted(self, '_log_cpms')
        X = check_array(X, dtype=np.int64)
        if hasattr(self, '_train_done'):
            memory = check_memory(self.memory)
            X = memory.cache(edger_tmm_logcpm_transform)(
                X, ref_sample=self.ref_sample_,
                prior_count=self.prior_count)[0]
        else:
            X = self._log_cpms
            self._train_done = True
        return super().transform(X)

    def inverse_transform(self, X):
        """
        Parameters
        ----------
        X : array-like, shape = (n_samples, n_features)
            The input data matrix.

        Returns
        -------
        X_r : array of shape (n_samples, n_original_features)
            `X` with columns of zeros inserted where features would have
            been removed by :meth:`transform`.
        """
        raise NotImplementedError("inverse_transform not implemented.")

    def _check_params(self, X, y):
        if not (self.k == "all" or 0 <= self.k <= X.shape[1]):
            raise ValueError(
                "k should be 0 <= k <= n_features; got %r."
                "Use k='all' to return all features."
                % self.k)

    def _get_support_mask(self):
        check_is_fitted(self, 'pvals_')
        mask = np.zeros_like(self.pvals_, dtype=bool)
        if self.k == 'all':
            mask = np.ones_like(self.pvals_, dtype=bool)
        elif self.k > 0:
            mask[np.argsort(self.pvals_, kind='mergesort')[:self.k]] = True
        return mask


class EdgeRFilterByExpr(BaseEstimator, SelectorMixin):
    """edgeR filterByExpr feature selector for RNA-seq data
    """
    def __init__(self):
        pass

    def fit(self, X, y):
        """
        Parameters
        ----------
        X : array-like, shape = (n_samples, n_features)
            The input counts data matrix.

        y : array-like, shape = (n_samples,)
            The target values (class labels in classification, real numbers in
            regression).

        Returns
        -------
        self : object
            Returns self.
        """
        X, y = check_X_y(X, y, dtype=np.int64)
        self._mask = np.array(r_edger_filterbyexpr_mask(X, y), dtype=bool)
        return self

    def _get_support_mask(self):
        check_is_fitted(self, '_mask')
        return self._mask


class LimmaVoom(BaseEstimator, SelectorMixin):
    """limma-voom feature selector and normalizer/transformer for RNA-seq data

    Parameters
    ----------
    k : int or "all" (default = "all")
        Number of top features to select. The "all" option bypasses selection,
        for use in a parameter search.

    lfc : float (default = 0)
        treat absolute log fold-change minimum threshold.  Default value
        of 0 gives eBayes results.

    robust : bool (default = True)
        limma treat/eBayes robust option

    prior_count : int (default = 1)
        Average count to add to each observation to avoid taking log of zero.
        Larger values for produce stronger moderation of the values for low
        counts and more shrinkage of the corresponding log fold changes.

    memory : None, str or object with the joblib.Memory interface \
        (default = None)
        Used for internal caching. By default, no caching is done.
        If a string is given, it is the path to the caching directory.

    Attributes
    ----------
    pvals_ : array, shape (n_features,)
        Feature raw p-values.

    padjs_ : array, shape (n_features,)
        Feature adjusted p-values.

    ref_sample_ : array, shape (n_features,)
        TMM normalization reference sample feature vector.
    """
    def __init__(self, k='all', lfc=0, robust=True, prior_count=1,
                 memory=None):
        self.k = k
        self.lfc = lfc
        self.robust = robust
        self.prior_count = prior_count
        self.memory = memory

    def fit(self, X, y):
        """
        Parameters
        ----------
        X : array-like, shape = (n_samples, n_features)
            The training input counts data matrix.

        y : array-like, shape = (n_samples,)
            The target values (class labels in classification, real numbers in
            regression).

        Returns
        -------
        self : object
            Returns self.
        """
        X, y = check_X_y(X, y, dtype=np.int64)
        self._check_params(X, y)
        memory = check_memory(self.memory)
        self.pvals_, self.padjs_, self._log_cpms, self.ref_sample_ = (
            memory.cache(limma_voom_feature_score)(
                X, y, lfc=self.lfc, robust=self.robust,
                prior_count=self.prior_count))
        return self

    def transform(self, X):
        """
        Parameters
        ----------
        X : array-like, shape = (n_samples, n_features)
            The input counts data matrix.

        Returns
        -------
        X_r : array of shape (n_samples, n_selected_features)
            edgeR TMM normalized log-CPM transformed input data matrix with
            only the selected features.
        """
        check_is_fitted(self, '_log_cpms')
        X = check_array(X, dtype=np.int64)
        if hasattr(self, '_train_done'):
            memory = check_memory(self.memory)
            X = memory.cache(edger_tmm_logcpm_transform)(
                X, ref_sample=self.ref_sample_,
                prior_count=self.prior_count)[0]
        else:
            X = self._log_cpms
            self._train_done = True
        return super().transform(X)

    def inverse_transform(self, X):
        """
        Parameters
        ----------
        X : array-like, shape = (n_samples, n_features)
            The input data matrix.

        Returns
        -------
        X_r : array of shape (n_samples, n_original_features)
            `X` with columns of zeros inserted where features would have
            been removed by :meth:`transform`.
        """
        raise NotImplementedError("inverse_transform not implemented.")

    def _check_params(self, X, y):
        if not (self.k == "all" or 0 <= self.k <= X.shape[1]):
            raise ValueError(
                "k should be 0 <= k <= n_features; got %r."
                "Use k='all' to return all features."
                % self.k)

    def _get_support_mask(self):
        check_is_fitted(self, 'pvals_')
        mask = np.zeros_like(self.pvals_, dtype=bool)
        if self.k == 'all':
            mask = np.ones_like(self.pvals_, dtype=bool)
        elif self.k > 0:
            mask[np.argsort(self.pvals_, kind='mergesort')[:self.k]] = True
        return mask


class LimmaScorerClassification(BaseScorer):
    """limma feature scorer for classification tasks.

    Parameters
    ----------
    robust : bool (default = False)
        limma eBayes robust option

    trend : bool (default = False)
        limma eBayes trend option

    Attributes
    ----------
    scores_ : array, shape (n_features,)
        Feature F values.

    pvalues_ : array, shape (n_features,)
        Feature adjusted p-values.
    """
    def __init__(self, robust=False, trend=False):
        self.robust = robust
        self.trend = trend

    def fit(self, X, y):
        """Run scorer on (X, y).

        Parameters
        ----------
        X : array_like, shape (n_samples, n_features)
            Feature matrix.

        y : array_like, shape (n_samples,)
            Target vector.

        Returns
        -------
        self : object
            Returns self.
        """
        self._check_params(X, y)
        f, pv = r_limma_feature_score(X, y, robust=self.robust,
                                      trend=self.trend)
        self.scores_ = np.array(f, dtype=float)
        self.pvalues_ = np.array(pv, dtype=float)
        return self

    def _check_params(self, X, y):
        pass


class ColumnSelector(BaseEstimator, SelectorMixin):
    """Manual column feature selector

    Parameters
    ----------
    cols : array-like (default = None)
        A list specifying the feature indices to be selected. For example,
        [1, 4, 5] to select the 2nd, 5th, and 6th feature columns.
        If None, returns all columns in the array.

    drop_axis : bool (default = False)
        Drops last axis if True and the only one column is selected. This
        is useful, e.g., when the ColumnSelector is used for selecting
        only one column and the resulting array should be fed to e.g.,
        a scikit-learn column selector. E.g., instead of returning an
        array with shape (n_samples, 1), drop_axis=True will return an
        aray with shape (n_samples,).
    """
    def __init__(self, cols=None, drop_axis=False):
        self.cols = cols
        self.drop_axis = drop_axis

    def fit(self, X, y):
        """
        Parameters
        ----------
        X : array-like, shape = (n_samples, n_features)
            The training input samples.

        y : array-like, shape = (n_samples)
            The target values (class labels in classification, real numbers in
            regression).

        Returns
        ---------
        self : object
            Returns self.
        """
        X, y = check_X_y(X, y, dtype=None)
        self._check_params(X, y)
        self._n_features = X.shape[1]
        return self

    def _check_params(self, X, y):
        if self.cols is not None:
            for col in self.cols:
                if not 0 <= col <= X.shape[1]:
                    raise ValueError(
                        "cols should be 0 <= col <= n_features; got %r."
                        "Use cols=None to return all features."
                        % col)

    def _get_support_mask(self):
        check_is_fitted(self, '_n_features')
        if self.cols is None:
            mask = np.ones(self._n_features, dtype=bool)
        else:
            mask = np.zeros(self._n_features, dtype=bool)
            mask[list(self.cols)] = True
        return mask


class CFS(BaseEstimator, SelectorMixin):
    """Feature selector using Correlation Feature Selection (CFS) algorithm

    Attributes
    ----------
    selected_idxs_ : array-like, 1d
        CFS selected feature indexes
    """
    def __init__(self):
        pass

    def fit(self, X, y):
        """
        Parameters
        ----------
        X : array-like, shape = (n_samples, n_features)
            The training input samples.

        y : array-like, shape = (n_samples)
            The target values (class labels in classification, real numbers in
            regression).

        Returns
        -------
        self : object
            Returns self.
        """
        X, y = check_X_y(X, y, dtype=None)
        self._n_features = X.shape[1]
        warnings.filterwarnings('ignore', category=RRuntimeWarning,
                                message="^Rjava\.init\.warning")
        self.selected_idxs_ = np.array(r_cfs_feature_idxs(X, y), dtype=int)
        warnings.filterwarnings('always', category=RRuntimeWarning)
        return self

    def _get_support_mask(self):
        check_is_fitted(self, 'selected_idxs_')
        mask = np.zeros(self._n_features, dtype=bool)
        mask[self.selected_idxs_] = True
        return mask


class FCBF(BaseEstimator, SelectorMixin):
    """Feature selector using Fast Correlation-Based Filter (FCBF) algorithm

    Parameters
    ----------
    k : int or "all" (default = "all")
        Number of top features to select. The "all" option bypasses selection,
        for use in a parameter search. If k is specified threshold is ignored.

    threshold : float (default = 0)
        The threshold value to use for feature selection. Features whose
        importance is greater or equal are kept while the others are
        discarded.

    memory : None, str or object with the joblib.Memory interface \
        (default = None)
        Used for internal caching. By default, no caching is done.
        If a string is given, it is the path to the caching directory.

    Attributes
    ----------
    selected_idxs_ : array-like, 1d
        FCBF selected feature indexes
    """
    def __init__(self, k='all', threshold=0, memory=None):
        self.k = k
        self.threshold = threshold
        self.memory = memory
        self.selected_idxs_ = np.array([], dtype=int)
        self.scores_ = np.array([], dtype=float)

    def fit(self, X, y):
        """
        Parameters
        ----------
        X : array-like, shape = (n_samples, n_features)
            The training input samples.

        y : array-like, shape = (n_samples)
            The target values (class labels in classification, real numbers in
            regression).

        Returns
        -------
        self : object
            Returns self.
        """
        X, y = check_X_y(X, y, dtype=None)
        self._check_params(X, y)
        memory = check_memory(self.memory)
        self._n_features = X.shape[1]
        if self.k == 'all' or self.k > 0:
            warnings.filterwarnings('ignore', category=RRuntimeWarning,
                                    message="^Rjava\.init\.warning")
            feature_idxs, scores = memory.cache(fcbf_feature_idxs)(
                X, y, threshold=self.threshold)
            warnings.filterwarnings('always', category=RRuntimeWarning)
            if self.k != 'all':
                feature_idxs = feature_idxs[
                    np.argsort(scores, kind='mergesort')[-self.k:]]
                scores = np.sort(scores, kind='mergesort')[-self.k:]
            self.selected_idxs_ = np.sort(feature_idxs, kind='mergesort')
            self.scores_ = scores[np.argsort(feature_idxs, kind='mergesort')]
        return self

    def _check_params(self, X, y):
        if not (self.k == "all" or 0 <= self.k <= X.shape[1]):
            raise ValueError(
                "k should be 0 <= k <= n_features; got %r."
                "Use k='all' to return all features."
                % self.k)

    def _get_support_mask(self):
        check_is_fitted(self, 'selected_idxs_')
        mask = np.zeros(self._n_features, dtype=bool)
        if self.k == 'all' or self.k > 0:
            mask[self.selected_idxs_] = True
        return mask


class ReliefF(BaseEstimator, SelectorMixin):
    """Feature selector using ReliefF algorithm

    Parameters
    ----------
    k : int or "all" (default = 10)
        Number of top features to select. The "all" option bypasses selection,
        for use in a parameter search.

    n_neighbors : int (default = 10)
        Number of neighbors for ReliefF algorithm

    sample_size : int (default = 5)
        Sample size for ReliefF algorithm

    memory : None, str or object with the joblib.Memory interface \
        (default = None)
        Used for internal caching. By default, no caching is done.
        If a string is given, it is the path to the caching directory.

    Attributes
    ----------
    scores_ : array-like, shape=(n_features,)
        Feature scores
    """
    def __init__(self, k=10, n_neighbors=10, sample_size=5, memory=None):
        self.k = k
        self.n_neighbors = n_neighbors
        self.sample_size = sample_size
        self.memory = memory

    def fit(self, X, y):
        """
        Parameters
        ----------
        X : array-like, shape = (n_samples, n_features)
            The training input samples.

        y : array-like, shape = (n_samples)
            The target values (class labels in classification, real numbers in
            regression).

        Returns
        -------
        self : object
            Returns self.
        """
        X, y = check_X_y(X, y, dtype=None)
        self._check_params(X, y)
        memory = check_memory(self.memory)
        warnings.filterwarnings('ignore', category=RRuntimeWarning,
                                message="^Rjava\.init\.warning")
        self.scores_ = memory.cache(relieff_feature_score)(X, y)
        warnings.filterwarnings('always', category=RRuntimeWarning)
        return self

    def _check_params(self, X, y):
        if not (self.k == "all" or 0 <= self.k <= X.shape[1]):
            raise ValueError(
                "k should be 0 <= k <= n_features; got %r."
                "Use k='all' to return all features."
                % self.k)

    def _get_support_mask(self):
        check_is_fitted(self, 'scores_')
        mask = np.zeros_like(self.scores_, dtype=bool)
        if self.k == 'all':
            mask = np.ones_like(self.scores_, dtype=bool)
        elif self.k > 0:
            mask[np.argsort(self.scores_, kind='mergesort')[-self.k:]] = True
        return mask
