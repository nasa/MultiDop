import numpy as np
import pyart

DEFAULT_KW = {
    'dir': './',
    'refl': 'DT',
    'vt': 'VT',
    'bgfile': None,
    'frprmn_out': 'frprmn_out.nc',
    'writeout': 'multidop_out.nc',
    'min_cba': 20.0,
    'calc_params': 'calc_params.dda',
    'anel': 1,
    'laplace': 1,
    'read_dataweights': 2,
    'read_cvg': 0,
    'max_dist': 10.0,
    'cutoff': 0.0,
    'UT': 0.0, 'VT': 0.0,
    'output_error': 1,
    'weak_height': -1,
    'upper_bc': 1,
    'itmax_frprmn': [200, 10],
    'itmax_dbrent': 200,
    'C1b': 1.0, 'C2b': 10.0, 'C3b': 0,
    'C4b': 1.0, 'C5b': 0.0, 'C8b': 0.0,
    'vary_weights': 0,
    'filter': ['none', '', ''],
    'cvg_opt_bg': [1, 1, 1], 'cvg_sub_bg': [0, 0, 0],
    'cvg_opt_fil': [0, 1, 1], 'cvg_sub_fil': [0, 0, 0],
    'cvg_bg': [0, 0, 0], 'cvg_fil': [0, 0, 0], 'sseq_trip': [1.0, 1.0, 0.0]}


class _Parameters(object):

    def __init__(self, params_dict):
        self.params = params_dict

    def gen_basic_string(self, var):
        return var + ':  ' + str(self.params[var]) + '\n'

    def gen_str_from_list(self, key):
        """
        Currently not clear whether this would be Python 3 compliant
        """
        outstr = key + ':  '
        for item in self.params[key]:
            outstr += str(item)+'  '
        outstr += '\n'
        return outstr

    def gen_radar_coverage(self, key, dual=True):
        nrads = len(self.params['files'])
        if not dual:
            outstr = ''
            for j in range(nrads):
                outstr += key + ':  ' + self.params['radar_names'][j] + \
                    '  ' + str(self.params[key][j]) + '\n'
            return outstr
        if nrads == 2:
            return key + ':  2  ' + self.params['radar_names'][0] + '  ' + \
                self.params['radar_names'][1] + '  ' + \
                str(self.params[key][0]) + '\n'
        elif nrads == 3:
            outstr = ''
            outstr += key + ':  2  ' + self.params['radar_names'][1] + \
                '  ' + self.params['radar_names'][2] + '  ' + \
                str(self.params[key][0]) + '\n'
            outstr += key + ':  2  ' + self.params['radar_names'][0] + \
                '  ' + self.params['radar_names'][2] + '  ' + \
                str(self.params[key][1]) + '\n'
            outstr += key + ':  2  ' + self.params['radar_names'][0] + '  ' + \
                self.params['radar_names'][1] + '  ' + \
                str(self.params[key][2]) + '\n'
            return outstr


class ParamFile(_Parameters):

    def __init__(self, params_dict, filename, ):
        if 'files' not in params_dict.keys() or \
                'radar_names' not in params_dict.keys():
            raise IOError('Must specify files and radar_names keys!')
        params_dict = check_kwargs(params_dict, DEFAULT_KW)
        _Parameters.__init__(self, params_dict)
        self.set_blockquotes()
        self.write_params_file(filename)

    def set_blockquotes(self):
        self.blockquotes = []
        self.blockquotes.append(
            """# Run time parameters for DDA.
# White space can be any sequence of spaces, tabs, and newlines.
# token: specifies type of parameter. Subsequent content depends on token.
# Parameters can be provided in any order.
# '#' to end of line is ignored.
""")
        self.blockquotes.append(
            """
# (Optional) Execution directory. Analysis process will set working directory
# to this when it starts. Auxiliary files hard coded into this program must be
# in this directory, or defaults will be used. Output files will appear here,
# too. Defaults to ".".
""")
        self.blockquotes.append(
            """
# Grid specification.
# axis: min step num_points
# max will be min + (num_points - 1) * step
# Assume, perhaps naively, that distances are measured in meters,
# angles in degrees.
""")
        self.blockquotes.append(
            """
# Geographic coordinates of grid origin, degrees
# Note that order is LONGITUDE LATITUDE, not vice versa (think Cartesian).
""")
        self.blockquotes.append(
            """
# OPAWS or OBAN input. Number of radars followed by that many file paths and
# radar name. Radar name identifies a radar elsewhere in parameter input.
""")
        self.blockquotes.append(
            """
# (Optional) Name of reflectivity variable. Defaults to "DT".
""")
        self.blockquotes.append(
            """
# (Optional) Name of velocity variable. Defaults to "VT".
""")
        self.blockquotes.append(
            """
# (Optional) Name of background file. If absent, analysis will not use sounding
# constraint.
""")
        self.blockquotes.append(
            """
# Names of output file
""")
        self.blockquotes.append(
            """
# Minimum beam crossing angle for optimal coverage
""")
        self.blockquotes.append(
            """
# Additional parameters for calculation
""")

    def write_params_file(self, filename):
        f = open(filename, 'w+')
        for i, bq in enumerate(self.blockquotes):
            f.write(bq)
            if i == 1:  # Path
                f.write('dir:  '+self.params['dir']+'\n')
            if i == 2:  # Grid specs
                for key in ['x', 'y', 'z']:
                    f.write(self.gen_str_from_list(key))
            if i == 3:
                f.write(self.gen_str_from_list('grid'))
            if i == 4:  # Radar files and names
                nrads = len(self.params['files'])
                f.write('opaws:  '+str(nrads)+'\n')
                for j in range(nrads):
                    f.write(self.params['files'][j] + ' ' +
                            self.params['radar_names'][j] + '\n')
            if i == 5:  # Reflectivity field name
                f.write(self.gen_basic_string('refl'))
            if i == 6:  # Reflectivity field name
                f.write(self.gen_basic_string('vt'))
            if i == 7:  # Background winds
                if self.params['bgfile'] is not None:
                    f.write(self.gen_basic_string('bgfile'))
            if i == 8:  # Output files
                f.write(self.gen_basic_string('frprmn_out'))
                f.write(self.gen_basic_string('writeout'))
            if i == 9:  # Beam crossing angle
                f.write(self.gen_basic_string('min_cba'))
            if i == 10:  # Calc Params file
                f.write(self.gen_basic_string('calc_params'))
        f.close()


class CalcParamFile(_Parameters):

    def __init__(self, params_dict, filename):
        _Parameters.__init__(self, params_dict)
        self.set_blockquotes()
        self.write_params_file(filename)

    def set_blockquotes(self):
        self.blockquotes = []
        self.blockquotes.append(
            """# For mass conservation constraint (CalcDiv):
# 0 => Boussinesq approximation
# 1 => anelastic approximation
""")
        self.blockquotes.append(
            """
# For smoothness constraint (CalcSmooth):
# 0 => first-order derivatives
# 1 => second-order derivatives
""")
        self.blockquotes.append(
            """
# For weights in data constraint:
# 0 => calculate and output to file
# 1 => read from file
# 2 => weight all observations equally
""")
        self.blockquotes.append(
            """
# For dual-Doppler coverage_bg mask
# Verification statistics computed within 2+ Doppler domain only.
# 0 => calculate
# 1 => read from file
# 2 => don't calculate
""")
        self.blockquotes.append(
            """
# Dual-Doppler domain criteria:
# There must be at least one observation from each radar within this distance
# of the analysis point (only used if read_cvg=0)
""")
        self.blockquotes.append(
            """
# Optional height below which to omit observations from analysis
# (e.g., data-denial experiments).
""")
        self.blockquotes.append(
            """
# Prescribed storm motion (spatiotemporally constant)
""")
        self.blockquotes.append(
            """
# If = 1, output verification stats after each iteration
# (see bottom of CalcCost)
""")
        self.blockquotes.append(
            """
# The index of the height equal to and below that
# the sounding constraint is weakened inside regions
# with greater than 10 dbz. Set to -1 to prevent
# implementation
""")
        self.blockquotes.append(
            """
# Impose upper BC w=0 if 1, ignore if -1
""")
        self.blockquotes.append(
            """
# Number of iterations in frprmn before and after the filtering round
""")
        self.blockquotes.append(
            """
# Number of iterations in dbrent
""")
        self.blockquotes.append(
            """
# Coefficient for data constraint
""")
        self.blockquotes.append(
            """
# Coefficient for mass conservation equation
""")
        self.blockquotes.append(
            """
# Coefficient for vorticity equation
""")
        self.blockquotes.append(
            """
# Coefficient for horizontal smoothing
""")
        self.blockquotes.append(
            """
# Coefficient for vertical smoothing
""")
        self.blockquotes.append(
            """
# Coefficient for sounding constraint
""")
        self.blockquotes.append(
            """
# 0 => all constant weights
# 1 => vary weights
""")
        self.blockquotes.append(
            """
# Define filter with ONE of the following forms.
# filter: none
# filter: filter_frequency Leise nstep
# filter: filter_frequency low-pass alpha
#
# Filter will be applied every filter_frequency iterations.
# Leise filter will perform nstep steps. See leise_filt.f
# Low pass filter will use alpha as smoothing parameter.
# Default is no filter.
""")
        self.blockquotes.append(
            """
# Coverage values for various combinations of radars.
# Each line should provide the type of coverage value, radar count,
# radar names, and the value, in the following form:
#
#   cvg_(""|opt|sub)_(bg|fil): integer radar1 radar2 ... boolean
#
# Radars are identified by the OPAWS/OBAN file name with grid data for that
# radar. This must be just the base name, not the full path.
#
# For example:
#
#   cvg_opt_bg: SR1 SR2 1
#
# says that if SR1 SR2
# both have data within max_dist meters of the point under consideration,
# and an optimal beam crossing angle, then the point will receive a coverage
# value of 1, i.e. point has coverage.
#
# "opt" means optimal beam crossing angle.
# "sub" means suboptimal beam crossing angle.
# "bg" means background coverage.
# "fil" means filter coverage.
# cvg_bg, cvg_fil, and sseq_trip do not require a radar count. (Beam crossing
# angle is meaningless with one radar, so there is no opt or sub)
#
# If this file is being used, coverage values must be provided for all
# combinations of radars.
""")
        self.blockquotes.append(
            """
# Background coverage arrays for optimal beam crossing angle for two radars
""")
        self.blockquotes.append(
            """
# Background coverage arrays for suboptimal beam crossing angle for two radars
""")
        self.blockquotes.append(
            """
# Filter coverage arrays for optimal beam crossing angle for two radars
""")
        self.blockquotes.append(
            """
# Filter coverage arrays for suboptimal beam crossing angle for two radars
""")
        self.blockquotes.append(
            """
# Background coverage arrays for single radar
""")
        self.blockquotes.append(
            """
# Filter coverage arrays for single radar
""")
        self.blockquotes.append(
            """
# Data weight multiplier when all 3 radar pairs have good beam crossing angles
""")

    def write_params_file(self, filename):
        f = open(filename, 'w+')
        for i, bq in enumerate(self.blockquotes):
            f.write(bq)
            if i == 0:
                f.write(self.gen_basic_string('anel'))
            if i == 1:
                f.write(self.gen_basic_string('laplace'))
            if i == 2:
                f.write(self.gen_basic_string('read_dataweights'))
            if i == 3:
                f.write(self.gen_basic_string('read_cvg'))
            if i == 4:
                f.write(self.gen_basic_string('max_dist'))
            if i == 5:
                f.write(self.gen_basic_string('cutoff'))
            if i == 6:
                f.write(self.gen_basic_string('UT'))
                f.write(self.gen_basic_string('VT'))
            if i == 7:
                f.write(self.gen_basic_string('output_error'))
            if i == 8:
                f.write(self.gen_basic_string('weak_height'))
            if i == 9:
                f.write(self.gen_basic_string('upper_bc'))
            if i == 10:
                f.write(self.gen_str_from_list('itmax_frprmn'))
            if i == 11:
                f.write(self.gen_basic_string('itmax_dbrent'))
            if i == 12:
                f.write(self.gen_basic_string('C1b'))
            if i == 13:
                f.write(self.gen_basic_string('C2b'))
            if i == 14:
                f.write(self.gen_basic_string('C3b'))
            if i == 15:
                f.write(self.gen_basic_string('C4b'))
            if i == 16:
                f.write(self.gen_basic_string('C5b'))
            if i == 17:
                f.write(self.gen_basic_string('C8b'))
            if i == 18:
                f.write(self.gen_basic_string('vary_weights'))
            if i == 19:
                f.write(self.gen_str_from_list('filter'))
            if i == 21:
                f.write(self.gen_radar_coverage('cvg_opt_bg'))
            if i == 22:
                f.write(self.gen_radar_coverage('cvg_sub_bg'))
            if i == 23:
                f.write(self.gen_radar_coverage('cvg_opt_fil'))
            if i == 24:
                f.write(self.gen_radar_coverage('cvg_sub_fil'))
            if i == 25:
                f.write(self.gen_radar_coverage('cvg_bg', dual=False))
            if i == 26:
                f.write(self.gen_radar_coverage('cvg_fil', dual=False))
            if i == 27:
                f.write(self.gen_radar_coverage('sseq_trip', dual=False))
        f.close()


def check_kwargs(kwargs, default_kw):
    """
    Check user-provided kwargs against defaults, and if some defaults aren't
    provided by user make sure they are provided to the DDA input file.

    Parameters
    ----------
    kwargs : dict
        User-specified keyword dictionary
    default_kw : dict
        Default keyword dictionary

    Returns
    -------
    kwargs : dict
        Dictionary of keywords for input to the DDA input parameter files
    """
    if 'grid' not in kwargs or 'x' not in kwargs or 'y' not in kwargs or \
            'z' not in kwargs:
        # Assumes matched grid in all files
        grid = pyart.io.read_grid(kwargs['files'][0])
    if 'grid' not in kwargs:
        kwargs['grid'] = [grid.origin_longitude['data'][0],
                          grid.origin_latitude['data'][0],
                          grid.origin_altitude['data'][0]]
    if 'x' not in kwargs:
        kwargs['x'] = [grid.x['data'].min(),
                       np.median(np.diff(grid.x['data'])),
                       len(grid.x['data'])]
    if 'y' not in kwargs:
        kwargs['y'] = [grid.y['data'].min(),
                       np.median(np.diff(grid.y['data'])),
                       len(grid.y['data'])]
    if 'z' not in kwargs:
        kwargs['z'] = [grid.z['data'].min(),
                       np.median(np.diff(grid.z['data'])),
                       len(grid.z['data'])]
    for key in default_kw:
        if key not in kwargs:
            kwargs[key] = default_kw[key]
    return kwargs
