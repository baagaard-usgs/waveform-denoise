#!/usr/bin/env python
#
# ======================================================================
#
#                           Brad T. Aagaard
#                        U.S. Geological Survey
#
# ======================================================================
#

import os
import logging
import numpy

DEFAULTS = """
[waveformsapp]
event = nc72282711
name = South Napa

[compare_figure]
height = 10.0
width = 7.5
margins = ((0.6, 0.5, 0.2), (0.5, 0.6, 0.7))


[files]
waveforms_acc = waveforms_acc.p
waveforms_vel = waveforms_vel.p
plots = plots-compare
"""

# ----------------------------------------------------------------------
def _config_get_list(list_string):
    """Convert list as string to list.

    :type list_string: list
    :param list_string: List as string.
    :returns: List of strings.
    """
    l = [f.strip() for f in list_string[1:-1].split(",")]
    return l

def _data_filename(params, pname, args=None, make_dir=False):
    """Construct relative path for file.
    
    :type params: ConfigParser
    :param params: Application parameters.
    :type pname: str
    :param pname: Name of parameter for file.
    :type args: tuple
    :param args: Tuple for arguments for substitution in name of parameter file.
    :type make_dir: bool
    :param make_dir: Create directory for file if True.
    """
    eventId = params.get("waveformsapp", "event")
    eventName = params.get("waveformsapp", "name").replace(" ","")
    eventDir = "%s-%s" % (eventId, eventName)

    if make_dir and not os.path.isdir(os.path.join("data", eventDir)):
        os.makedirs(os.path.join("data", eventDir))

    filename = params.get("files", pname) if args is None else params.get("files", pname) % args
    return os.path.join("data", eventDir, filename)
    

# ----------------------------------------------------------------------
class WaveformData(object):
    """Object for retrieving and processing waveform data.
    """

    def __init__(self, params, show_progress=True):
        """Constructor.

        :type params: ConfigParser
        :param params: Parameters for application.
        :type show_progress: bool
        :param show_progress: Show progress during execution.
        """
        self.params = params
        self.showProgress = show_progress

        self.accSM = None
        self.velBB = None
        return

    def load_processed_waveforms(self):
        """Load processed data.
        """
        import cPickle
        
        if self.showProgress:
            print("Loading processed waveforms...")            

        if self.accSM is None:
            with open(_data_filename(self.params, "waveforms_acc"), "r") as fin:
                self.accSM = cPickle.Unpickler(fin).load()
        if self.velBB is None:
            with open(_data_filename(self.params, "waveforms_vel"), "r") as fin:
                self.velBB = cPickle.Unpickler(fin).load()

        return
    
# ----------------------------------------------------------------------
class ComparisonFigure(object):
    """
    Rows are comprison and residual for original and denoised.
    Columns are E, N, and Z components
    """

    ROWS = ["Bandpassed", "Bandpassed Residual", "Denoised", "Denoised Residual"]
    COLS = ["E", "N", "Z"]
    
    def __init__(self, params, showProgress):
        """
        """
        self.params = params
        self.showProgress = showProgress
        return


    def plot(self, data):
        """
        """
        import sys
        from ast import literal_eval
        import obspyutils.subset
        import obspyutils.baseline
        
        if self.showProgress:
            sys.stdout.write("Plotting comparison figures...")

        from basemap.Figure import Figure
        
        self.figure = Figure()
        w = self.params.getfloat("compare_figure", "width")
        h = self.params.getfloat("compare_figure", "height")
        margins = literal_eval(self.params.get("compare_figure", "margins"))
        self.figure.open(w, h, margins)

        self._setupSubplots()
        self.figure.figure.canvas.draw()

        accSM = obspyutils.subset.streamByStation(data.accSM["bandpass"])
        velBB = obspyutils.subset.streamByStation(data.velBB["orig"])
        for st in accSM.keys():
            if not st in velBB:
                continue
            for tr in accSM[st]:
                samplingRate = tr.stats.sampling_rate
                startTime = tr.stats.starttime
                trBB = velBB[st].select(component=tr.stats.channel[-1])
                trBB.interpolate(sampling_rate=samplingRate, method="weighted_average_slopes", starttime=startTime)
                trBB.trim(starttime=tr.stats.starttime, endtime=tr.stats.endtime)

        accSM = obspyutils.subset.streamByStation(data.accSM["data"])
        accSMOrig = obspyutils.subset.streamByStation(data.accSM["bandpass"])
        velBB = obspyutils.subset.streamByStation(data.velBB["data"])
        velBBOrig = obspyutils.subset.streamByStation(data.velBB["orig"])

        numStations = len(accSM.keys())
        for ist,st in enumerate(accSM.keys()):
            if not st in velBB:
                continue
            
            info = "%s.%s" % (accSM[st].traces[0].stats.network, accSM[st].traces[0].stats.station)
            self.figure.figure.suptitle(info, fontweight='bold')

            for component in ["E", "N", "Z"]:
                trBB = velBBOrig[st].select(component=component)[0]

                # Original strong-motion versus broadband
                trSM = accSMOrig[st].select(component=component)[0]
                t = trSM.times()#reftime=originTime)
                velSM, dispSM = obspyutils.baseline.integrate_acc(trSM)

                self._updatePlot("Bandpassed", component, t, velSM.data, trBB.data)
                self._updateResidual("Bandpassed", component, t, velSM.data-trBB.data)

                # Denoised strong-motion versus broadband
                trSM = accSM[st].select(component=component)[0]
                t = trSM.times()#reftime=originTime)
                velSM, dispSM = obspyutils.baseline.integrate_acc(trSM)
                trBB.trim(starttime=velSM.stats.starttime, endtime=velSM.stats.endtime, pad=True, fill_value=0.0)
                
                self._updatePlot("Denoised", component, t, velSM.data, trBB.data)
                self._updateResidual("Denoised", component, t, velSM.data-trBB.data)
                

            plotsDir = os.path.join(_data_filename(self.params, "plots"))
            if not os.path.isdir(plotsDir):
                os.makedirs(plotsDir)
            self.figure.figure.savefig(os.path.join(plotsDir, info+".png"))

            if self.showProgress:
                sys.stdout.write("\rPlotting comparison figures...%d%%" % (((ist+1)*100)/numStations))
                sys.stdout.flush()
        if self.showProgress:
            sys.stdout.write("\n")
            
        return

    def _setupSubplots(self):
        """
        """
        nrows = len(ComparisonFigure.ROWS)
        ncols = len(ComparisonFigure.COLS)
        
        self.axes = {}
        for irow,row in enumerate(ComparisonFigure.ROWS):
            for icol,col in enumerate(ComparisonFigure.COLS):
                ax = self.figure.axes(nrows, ncols, irow+1, icol+1)
                line, = ax.plot([], [], 'k-', lw=0.5)
                line2, = ax.plot([], [], 'r-', lw=0.5)
                ax.autoscale(enable=True, axis="both", tight=True)
                if irow == 0:
                    ax.set_title(col)
                if irow == nrows-1:
                    ax.set_xlabel("Time (s)")
                if icol == 0:
                    ax.set_ylabel("Velocity (m/s)")
                    pos = ax.get_position()
                    ax.text(pos.xmin, pos.ymax+0.02, row, fontweight='bold', transform=self.figure.figure.transFigure, ha="left")
                self.axes["%s_%s" % (row,col)] = (ax,line,line2)

        return

    def _updatePlot(self, row, component, t, dSM, dBB):
        ax, lineBB, lineSM = self.axes[row+"_"+component]
        lineSM.set_xdata(t)
        lineSM.set_ydata(dSM)
        lineBB.set_xdata(t)
        lineBB.set_ydata(dBB)
        ax.relim()
        ax.autoscale_view()
        ax.draw_artist(lineSM)
        ax.draw_artist(lineBB)
        return
    
    def _updateResidual(self, row, component, t, residual):
        ax, line, lineR = self.axes[row+" Residual_"+component]
        lineR.set_xdata(t)
        lineR.set_ydata(residual)
        ax.relim()
        ax.autoscale_view()
        ax.draw_artist(lineR)
        return
    
# ----------------------------------------------------------------------
class ComparisonApp(object):
    """
    Plot comparison of denoised strong-motion and broadband data velocity waveforms.
    """
    
    def __init__(self, show_progress=True):
        """Constructor.

        :type show_progress: bool
        :param show_progress: Show progress during execution.
        """
        self.showProgress = show_progress
        return

    def initialize(self, config_filenames):
        """Set parameters from config file and DEFAULTS.

        :type config_filename: str
        :param config_filename: Name of configuration (INI) file with parameters.
        """
        import ConfigParser
        import io
        config = ConfigParser.SafeConfigParser()
        config.readfp(io.BytesIO(DEFAULTS))
        for filename in config_filenames.split(","):
            if self.showProgress:
                print("Fetching parameters from %s..." % filename)
            config.read(filename)

        self.params = config
        return
    
    def show_parameters(self):
        """Write parameters to stdout.
        """
        import sys
        self.params.write(sys.stdout)
        return

    def event_name(self):
        """Get event name from event id and earthquake name/location.

        :return: Event name
        """
        eventid = self.params.get("waveformsapp", "event")
        eventname = self.params.get("waveformsapp", "name").replace(" ","").lower()
        return "%s-%s" % (eventid, eventname)
        

# ======================================================================
if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser()
    parser.add_argument("--config", action="store", dest="config", required=True)
    parser.add_argument("--show-parameters", action="store_true", dest="show_parameters")
    parser.add_argument("--plot-comparison", action="store_true", dest="plot_comparison")
    parser.add_argument("--compute-metrics", action="store_true", dest="compute_metrics")
    parser.add_argument("--all", action="store_true", dest="all")
    parser.add_argument("--quiet", action="store_false", dest="show_progress", default=True)
    parser.add_argument("--debug", action="store_true", dest="debug")
    args = parser.parse_args()

    app = ComparisonApp(args.show_progress)
    app.initialize(args.config)

    logLevel = logging.DEBUG if args.debug else logging.INFO
    logging.basicConfig(level=logLevel, filename=app.event_name()+".log")

    data = WaveformData(app.params, args.show_progress)
    
    if args.show_parameters or args.all:
        app.show_parameters()
    
    if args.compute_metrics or args.all:
        data.compute_metrics()

    if args.plot_comparison or args.all:
        plotter = ComparisonFigure(app.params, args.show_progress)
        data.load_processed_waveforms()
        plotter.plot(data)
        
# End of file
