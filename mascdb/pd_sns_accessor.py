# -----------------------------------------------------------------------------.
# Copyright (c) 2021-2025 MASCDB developers
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the MIT License.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
#
# You should have received a copy of the MIT License
# along with this program.  If not, see <https://opensource.org/license/mit/>.
# -----------------------------------------------------------------------------.
"""Pandas DataFrame Seaborn Accessor."""
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns


def _get_sns_fun_names():
    sns_fun_names = [
        "boxplot",
        "violinplot",
        "boxplot",
        "boxenplot",
        "swarmplot",
        "stripplot",
        "pointplot",
        "lmplot",
        "pairplot",
        "scatterplot",
        "relplot",
        "lineplot",
        "displot",
        "catplot",
        "barplot",
        "histplot",
        "jointplot",
        "kdeplot",
    ]
    return sns_fun_names


def _sns_fun_factory(sns_fun_name):
    fun = getattr(sns, sns_fun_name)

    def method_fun(self, **kwargs):
        # Remove 'data' argument if specified
        _ = kwargs.pop("data", None)
        # Plot
        return fun(data=self._obj, **kwargs)

    return method_fun


@pd.api.extensions.register_dataframe_accessor("sns")
class SeabornAccessor:
    """Pandas DataFrame accessor for Seaborn plotting functionality.

    This accessor provides convenient access to Seaborn plotting functions directly
    from pandas DataFrames using the `.sns` attribute. It dynamically adds all common
    Seaborn plotting methods and provides additional custom visualization methods.

    Examples
    --------
    >>> import pandas as pd
    >>> df = pd.DataFrame({"x": [1, 2, 3], "y": [4, 5, 6]})
    >>> df.sns.scatterplot(x="x", y="y")
    >>> df.sns.corrplot()

    Notes
    -----
    Available Seaborn methods include: boxplot, violinplot, boxenplot, swarmplot,
    stripplot, pointplot, lmplot, pairplot, scatterplot, relplot, lineplot,
    displot, catplot, barplot, histplot, jointplot, and kdeplot.

    """

    def __init__(self, pandas_obj):
        self._obj = pandas_obj
        # Add all methods dynamically
        sns_fun_names = _get_sns_fun_names()
        for sns_fun_name in sns_fun_names:
            setattr(type(self), sns_fun_name, _sns_fun_factory(sns_fun_name))

    def corrplot(self, vars=None, vmin=-0.3, vmax=0.3, center=0, cbar_kws=None, linewidths=0.5):
        """Create a correlation matrix heatmap with lower triangle display.

        Computes the correlation matrix of the DataFrame and visualizes it as a heatmap
        with a mask for the upper triangle, showing only the lower triangle of correlations.

        Parameters
        ----------
        vars : list of str, optional
            List of column names to include in the correlation matrix.
            If None, uses all numeric columns. Default is None.
        vmin : float, optional
            Minimum value for colormap normalization. Default is -0.3.
        vmax : float, optional
            Maximum value for colormap normalization. Default is 0.3.
        center : float, optional
            Value at which to center the colormap. Default is 0.
        cbar_kws : dict, optional
            Keyword arguments for the colorbar. Default is {"shrink": 0.5}.
        linewidths : float, optional
            Width of lines separating cells in the heatmap. Default is 0.5.

        Returns
        -------
        matplotlib.figure.Figure
            The matplotlib figure object containing the correlation heatmap.

        Examples
        --------
        >>> df.sns.corrplot()
        >>> df.sns.corrplot(vars=["col1", "col2", "col3"], vmin=-1, vmax=1)

        References
        ----------
        https://seaborn.pydata.org/examples/many_pairwise_correlations.html

        """
        # https://seaborn.pydata.org/examples/many_pairwise_correlations.html
        if cbar_kws is None:
            cbar_kws = {"shrink": 0.5}
        df = self._obj.loc[:, vars] if vars is not None else self._obj
        # Compute the correlation matrix
        corr = df.corr()

        # Generate a mask for the upper triangle
        mask = np.triu(np.ones_like(corr, dtype=bool))

        # Set up the matplotlib figure
        f, ax = plt.subplots(figsize=(11, 9))

        # Generate a custom diverging colormap
        cmap = sns.diverging_palette(230, 20, as_cmap=True)

        # Draw the heatmap with the mask and correct aspect ratio
        sns.heatmap(
            corr,
            mask=mask,
            cmap=cmap,
            vmin=vmin,
            vmax=vmax,
            center=center,
            cbar_kws=cbar_kws,
            square=True,
            linewidths=linewidths,
        )
        return f

    def kde_marginals(
        self,
        x,
        y,
        xlim=None,
        ylim=None,
        space=0,
        thresh=0,
        levels=100,
        cmap="rocket",
        hist_color="#03051A",
        hist_alpha=1,
        hist_bins=25,
    ):
        """Create a bivariate KDE plot with marginal histograms.

        Produces a joint plot showing a smooth bivariate kernel density estimate (KDE)
        in the center with marginal histograms along the x and y axes.

        Parameters
        ----------
        x : str
            Name of the column to plot on the x-axis.
        y : str
            Name of the column to plot on the y-axis.
        xlim : tuple of float, optional
            Limits for the x-axis as (min, max). Default is None.
        ylim : tuple of float, optional
            Limits for the y-axis as (min, max). Default is None.
        space : float, optional
            Space between the joint and marginal axes. Default is 0.
        thresh : float, optional
            Threshold for the KDE contours. Values below this are not drawn. Default is 0.
        levels : int, optional
            Number of contour levels for the KDE plot. Default is 100.
        cmap : str, optional
            Colormap name for the KDE plot. Default is "rocket".
        hist_color : str, optional
            Color for the marginal histograms. Default is "#03051A".
        hist_alpha : float, optional
            Alpha (transparency) for the marginal histograms. Default is 1.
        hist_bins : int, optional
            Number of bins for the marginal histograms. Default is 25.

        Returns
        -------
        seaborn.JointGrid
            The JointGrid object containing the bivariate and marginal plots.

        Examples
        --------
        >>> df.sns.kde_marginals(x="col1", y="col2")
        >>> df.sns.kde_marginals(x="col1", y="col2", cmap="viridis", hist_bins=50)

        References
        ----------
        https://seaborn.pydata.org/examples/smooth_bivariate_kde.html

        """
        # https://seaborn.pydata.org/examples/smooth_bivariate_kde.html
        df = self._obj
        g = sns.JointGrid(data=df, x=x, y=y, space=space)
        g.plot_joint(sns.kdeplot, fill=True, clip=(xlim, ylim), thresh=thresh, levels=levels, cmap=cmap)
        g.plot_marginals(sns.histplot, bins=hist_bins, color=hist_color, alpha=hist_alpha)
        return g

    def kde_ridgeplot(self, x, group, pal=None, bw_adjust=0.5, height=0.5, aspect=15, hspace=-0.25, linewidth=2):
        """Create a ridge plot (joyplot) showing KDE distributions for groups.

        Produces overlapping kernel density estimates for different groups, creating
        a "ridge" or "joyplot" visualization useful for comparing distributions across
        multiple categories.

        Parameters
        ----------
        x : str
            Name of the column containing the values to plot.
        group : str
            Name of the column containing the grouping variable.
        pal : list or palette, optional
            Color palette for the groups. If None, uses a cubehelix palette.
            Default is None.
        bw_adjust : float, optional
            Bandwidth adjustment factor for the KDE. Higher values produce smoother curves.
            Default is 0.5.
        height : float, optional
            Height of each facet in inches. Default is 0.5.
        aspect : float, optional
            Aspect ratio of each facet (width/height). Default is 15.
        hspace : float, optional
            Space between facets. Negative values create overlap. Default is -0.25.
        linewidth : float, optional
            Width of the KDE lines. Default is 2.

        Returns
        -------
        seaborn.FacetGrid
            The FacetGrid object containing the ridge plot.

        Examples
        --------
        >>> df.sns.kde_ridgeplot(x="value", group="category")
        >>> df.sns.kde_ridgeplot(x="value", group="category", bw_adjust=1.0, hspace=-0.5)

        References
        ----------
        https://seaborn.pydata.org/examples/kde_ridgeplot.html

        """
        # https://seaborn.pydata.org/examples/kde_ridgeplot.html
        # sns.set_theme(style="white", rc={"axes.facecolor": (0, 0, 0, 0)})

        # Reorder dataframe by group
        df = self._obj
        df = df[[x, group]]
        # m = df[group].map(ord)
        # df["x"] += m

        # Initialize the FacetGrid object
        if pal is None:
            pal = sns.cubehelix_palette(10, rot=-0.25, light=0.7)

        g = sns.FacetGrid(df, row=group, hue=group, aspect=aspect, height=height, palette=pal)

        # Draw the densities in a few steps
        g.map(
            sns.kdeplot,
            x,
            bw_adjust=bw_adjust,
            linewidth=linewidth,
            clip_on=False,
            fill=True,
            alpha=1,
        )

        g.map(
            sns.kdeplot,
            x,
            bw_adjust=bw_adjust,
            lw=2,
            clip_on=False,
            color="w",
        )

        # Passing color=None to refline() uses the hue mapping
        g.refline(y=0, linewidth=2, linestyle="-", color=None, clip_on=False)

        # # Define and use a simple function to label the plot in axes coordinates
        # def label(x, color, label):
        #     ax = plt.gca()
        #     ax.text(0, 0.2, label, fontweight="bold", color=color, ha="left", va="center", transform=ax.transAxes)

        # g.map(label, x)

        def label(x, y, **kwargs):  # noqa: ARG001
            ax = plt.gca()
            #  kwargs contains things like color, label, etc.
            ax.text(
                0.05,
                0.2,
                kwargs.get("label", ""),
                fontweight="bold",
                color=kwargs.get("color"),
                ha="left",
                va="center",
                transform=ax.transAxes,
            )

        g.map(label)

        # Set the subplots to overlap
        g.figure.subplots_adjust(hspace=hspace)

        # Remove axes details that don't play well with overlap
        g.set_titles("")
        g.set(yticks=[], ylabel="")
        g.despine(bottom=True, left=True)

        # Return the object
        return g
