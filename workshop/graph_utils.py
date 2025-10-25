"""
Live plotting utilities for IMD streaming workshop.

This module provides a simple interface for creating live-updating plots
during MD simulation streaming with MDAnalysis and IMDv3.
"""

import matplotlib.pyplot as plt
from IPython.display import display, clear_output
import numpy as np


def live_plot(title="Live Analysis", xlabel="Time (ps)", ylabel="Value", 
              figsize=(10, 6), update_interval=1):
    """
    Create a live-updating plot for streaming MD analysis.
    
    Parameters
    ----------
    title : str
        Plot title
    xlabel : str
        X-axis label
    ylabel : str  
        Y-axis label
    figsize : tuple
        Figure size (width, height) in inches
    update_interval : int
        Update plot every N frames (1 = every frame, 10 = every 10th frame)
    
    Returns
    -------
    dict
        Dictionary with 'fig', 'ax', 'line', 'times', 'values' and 'update' function
        
    Examples
    --------
    >>> plot = live_plot("Distance", ylabel="Distance (Å)")
    >>> for ts in u.trajectory:
    >>>     distance = calculate_distance()
    >>>     plot['update'](ts.time, distance)
    """
    # Setup figure
    plt.ion()
    fig, ax = plt.subplots(figsize=figsize)
    line, = ax.plot([], [], 'b-', linewidth=2, label=ylabel)
    
    ax.set_xlabel(xlabel, fontsize=12)
    ax.set_ylabel(ylabel, fontsize=12)
    ax.set_title(title, fontsize=14, fontweight='bold')
    ax.legend(loc='best')
    ax.grid(True, alpha=0.3, linestyle='--')
    
    # Data storage
    times = []
    values = []
    frame_counter = [0]  # Use list to allow modification in nested function
    
    def update(time, value):
        """Update the plot with new data point."""
        frame_counter[0] += 1
        times.append(time)
        values.append(value)
        
        # Only update display every N frames for performance
        if frame_counter[0] % update_interval == 0:
            line.set_data(times, values)
            ax.relim()
            ax.autoscale_view()
            clear_output(wait=True)
            display(fig)
    
    # Return everything in a dict for easy access
    return {
        'fig': fig,
        'ax': ax,
        'line': line,
        'times': times,
        'values': values,
        'update': update,
        'close': lambda: plt.ioff()
    }


def live_plot_multi(title="Live Analysis", xlabel="Time (ps)", 
                    ylabels=None, figsize=(10, 6), update_interval=1):
    """
    Create a live-updating plot with multiple y-values for streaming MD analysis.
    
    Parameters
    ----------
    title : str
        Plot title
    xlabel : str
        X-axis label
    ylabels : list of str
        Labels for each y-value series
    figsize : tuple
        Figure size (width, height) in inches
    update_interval : int
        Update plot every N frames
    
    Returns
    -------
    dict
        Dictionary with 'fig', 'ax', 'lines', 'times', 'values_list' and 'update' function
        
    Examples
    --------
    >>> plot = live_plot_multi("Distances", ylabels=["N-ter", "C-ter"])
    >>> for ts in u.trajectory:
    >>>     dist1 = calculate_distance1()
    >>>     dist2 = calculate_distance2()
    >>>     plot['update'](ts.time, [dist1, dist2])
    """
    if ylabels is None:
        ylabels = ["Value"]
    
    # Setup figure
    plt.ion()
    fig, ax = plt.subplots(figsize=figsize)
    
    colors = plt.cm.tab10(np.linspace(0, 1, len(ylabels)))
    lines = []
    for i, (label, color) in enumerate(zip(ylabels, colors)):
        line, = ax.plot([], [], linewidth=2, label=label, color=color)
        lines.append(line)
    
    ax.set_xlabel(xlabel, fontsize=12)
    ax.set_ylabel("Value", fontsize=12)
    ax.set_title(title, fontsize=14, fontweight='bold')
    ax.legend(loc='best')
    ax.grid(True, alpha=0.3, linestyle='--')
    
    # Data storage
    times = []
    values_list = [[] for _ in ylabels]
    frame_counter = [0]
    
    def update(time, values):
        """Update the plot with new data points.
        
        Parameters
        ----------
        time : float
            Current time point
        values : list
            List of values corresponding to each y-series
        """
        frame_counter[0] += 1
        times.append(time)
        
        for i, val in enumerate(values):
            values_list[i].append(val)
        
        # Only update display every N frames for performance
        if frame_counter[0] % update_interval == 0:
            for i, line in enumerate(lines):
                line.set_data(times, values_list[i])
            ax.relim()
            ax.autoscale_view()
            clear_output(wait=True)
            display(fig)
    
    return {
        'fig': fig,
        'ax': ax,
        'lines': lines,
        'times': times,
        'values_list': values_list,
        'update': update,
        'close': lambda: plt.ioff()
    }
