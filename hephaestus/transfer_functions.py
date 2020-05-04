"""
Python module for computing thermocouple emf values from temperatures.
This module just contains the generic thermocouple class and helper
functions.
"""

__author__ = "lmbacelar@gmail.com"
__copyright__ = "public domain"

import numpy as np

# scipy.optimize will be imported when needed.
optimize = None


def ensure_import_optimize():
    global optimize
    if optimize == None:
        try:
            import scipy.optimize as optimize
        except ImportError:
            raise ImportError(
                "Inverse lookup requires scipy.optimize module. Please install SciPy."
            )


class Polynomial_Gaussian(object):
    """\
    Piecewise mathematical function of polynomials plus gaussian, used for
    thermocouple reference.
    
    Main methods:
     func(T)          # compute the function
     func.__call__(T) # synonym for func(T)
     func.inverse(F)  # perform inverse lookup
    
    The raw function parameters are stored in .table. The structure of .table
    is a list of tuples giving the different segments of the piecewise function,
    formatted as:
        (minimum T, maximum T, polynomial coefs array, exponential coefs list)
    The polynomial coefs array is in the order of np.polyval(), i.e., starting
    with the highest power and ending with zeroth power (offset).
    Exponential coefs are used as  ec[0] * np.exp(ec[1] * (T - ec[2])**2), or
    may be None in which case only the polynomial is used.
    
    The appropriate temperature and voltage units to use when calling these
    functions are assumed to be degrees Celsius and milivolts as defined in
    NIST, ASTM and OMEGA reference tables.
    
    .source and .calibration are strings containing information about where the
    function data comes from, and how it is calibrated.
    """

    def __init__(self, table=None, invtable=None, source="", calibration=""):
        self.table = table
        self.invtable = invtable
        self.source = source
        self.calibration = calibration

        # check table
        lastmax = table[0][0]
        for tmin, tmax, _, _ in table:
            if not tmin <= tmax:
                raise ValueError("Temperature limits must be in ascending order.")
            if tmin != lastmax:
                raise ValueError("Pieces' limits must be contiguous.")
            lastmax = tmax

    @property
    def minT(self):
        return self.table[0][0]

    @property
    def maxT(self):
        return self.table[-1][1]

    @property
    def minV(self):
        return self.invtable[0][0]

    @property
    def maxV(self):
        return self.invtable[-1][1]

    def __repr__(self):
        return (
            "<piecewise polynomial+gaussian, domain %g to %g in ºC, output in mV; %s calibrated, from %s>"
            % (self.minT, self.maxT, self.calibration, self.source,)
        )

    def __call__(self, T, derivative=0, out_of_range="raise"):
        """\
        Calculate reference function at given temperature.

        Parameters
        ----------
        T : array_like
            Temperature or array of temperatures.
        derivative: integer
            Use this parameter to evaluate the functional derivative of the emf
            function at a given temperature. Default is derivative=0 (no derivative).
        out_of_range: string, optional
            Determines behaviour for out of range temperatures.
            "raise": raises an ValueError exception. (default)
            "nan":   values replaced by nans.
            "extrapolate": extrapolates from closest range. Do not trust this!
        
        Returns
        -------
        emf : array_like
            computed emf function
        """
        if out_of_range not in ["raise", "nan", "extrapolate"]:
            raise ValueError("invalid out_of_range parameter", out_of_range)

        T = np.array(T, copy=False, order="A")
        emf_choices = [None]

        # We go through the table, determining the selector which is used
        # to choose which piece of the piecewise function to use.
        # selector = 0 where T is underrange,
        # selector = 1 where T is in first range,
        #  ...
        # selector = N where T is in last (Nth) range,
        # selector = N+1 where T is overrange.
        tmin = self.minT
        selector = (T >= tmin) * 1
        for tmin, tmax, coefs, ec in self.table:
            selector += T > tmax
            # Here we go ahead and compute emf values using all ranges.
            #   this is simple but perhaps a bit inefficient.
            emf = np.polyval(np.polyder(coefs, derivative), T)

            if ec:
                # Type K thermocouple has this annoying exponential addition term,
                # corresponding to a little bump at 127 Celsius.
                dT = T - ec[2]
                gauss = ec[0] * np.exp(ec[1] * dT ** 2)
                if derivative == 0:
                    emf += gauss
                elif derivative == 1:
                    emf += 2.0 * ec[1] * gauss * dT
                elif derivative == 2:
                    emf += 2.0 * ec[1] * gauss * (2.0 * ec[1] * dT ** 2 + 1.0)
                elif derivative == 3:
                    emf += (
                        4.0 * ec[1] * ec[1] * gauss * dT * (2.0 * ec[1] * dT ** 2 + 3.0)
                    )
                else:
                    raise ValueError(
                        "sorry, derivatives > 3 not supported for this type."
                    )
            emf_choices.append(emf)
        emf_choices.append(None)

        if out_of_range == "nan":
            emf_choices[0] = T * np.nan
            emf_choices[-1] = emf_choices[0]
        else:
            emf_choices[0] = emf_choices[1]
            emf_choices[-1] = emf_choices[-2]

        if out_of_range == "raise":
            unders = selector <= 0
            overs = selector > len(self.table)
            if np.any(unders) or np.any(overs):
                u_temps = np.extract(unders, T)
                o_temps = np.extract(overs, T)
                if u_temps.size == 0:
                    u_temps = None
                if o_temps.size == 0:
                    o_temps = None
                msg = "Temperatures (ºC) under or over range:"
                raise ValueError(msg, u_temps, o_temps)

        return np.choose(selector, emf_choices)

    def refinv(self, V):
        """\
        Calculate temperature at given voltage using reference 
        inverse polynomial function.

        Parameters
        ----------
        V : array_like
            Voltage or array of voltages.
        
        Returns
        -------
        t : array_like
            computed t function
        """
        V = np.array(V, copy=False, order="A")
        t_choices = [None]

        vmin = self.minV
        selector = (V >= vmin) * 1
        for vmin, vmax, coefs, _ in self.invtable:
            selector += V > vmax
            t = np.polyval(coefs, V)
            t_choices.append(t)
        t_choices.append(None)
        t_choices[0] = t_choices[1]
        t_choices[-1] = t_choices[-2]
        return np.choose(selector, t_choices)

    def inverse(self, V, Tstart=None, Vtol=1e-6):
        """
        Find the temperature corresponding to a given voltage, via zero-finding.
        
        Parameters
        ----------
        V: float
            Measured voltage (in milivolts) goes here.
        Tstart: float
            Suggested starting temperature for search. Defaults to reference
            inverse function or midpoint of range.
        Vtol: float
            Desired absolute tolerance of voltage value.
        
        Returns
        -------
        T: float
            Temperature T, such that func(T) = V
            If the solution does not converge within |func(T) - V| > Vtol,
            an exception is raised.

        Note on implementation
        ----------------------
        First checks if func(Tstart) is close enough to V;
        If this fails, try to use scipy.optimize.newton;
        Failing that, use scipy.optimize.brentq.
        
        This function requires scipy to be installed when using scipy.optimize.
        It will attemp to import it upon the first usage.
        """
        V = float(V)
        if Tstart == None:
            if self.invtable == None:
                Tstart = 0.5 * (self.minT + self.maxT)
            else:
                Tstart = self.refinv(V)

        if abs(self(Tstart, out_of_range="extrapolate") - V) <= Vtol:
            return Tstart

        ensure_import_optimize()
        fun0 = lambda T: self(T, out_of_range="extrapolate") - V
        fun1 = lambda T: self(T, derivative=1, out_of_range="extrapolate")
        fun2 = lambda T: self(T, derivative=2, out_of_range="extrapolate")
        try:
            T = optimize.newton(fun0, Tstart, fprime=fun1, fprime2=fun2, tol=Vtol)
            if abs(self(T, out_of_range="extrapolate") - V) > Vtol:
                raise ValueError
        except:
            try:
                T = optimize.brentq(fun0, self.minT, self.maxT)
            except ValueError as e:
                if e.args == ("f(a) and f(b) must have different signs",):
                    raise ValueError("Voltage not within in allowed range.")
                else:
                    raise
            if not abs(self(T, out_of_range="extrapolate") - V) <= Vtol:
                raise ValueError("Did not converge within tolerance.")

        return T


class Thermocouple(object):
    """
    Thermocouple helper object. This object provides practical
    methods for converting between temperatures and measured voltages:
    
    * ``.emf(T)`` returns voltage from known temperature.
    * ``.t(V)`` returns temperature from known voltage.
    
    Units according to reference tables - milivolt and degree Celsius

    In each case it is possible (and desirable) to pass in the reference
    junction temperature by the keyword argument Tref.
    """

    def __init__(self, func, ttype=""):
        """
        func is the object that contains the actual function information, and has
        methods __call__, inverse, and attributes .minT, .maxT.
        """
        self.func = func
        self.type = ttype

    def __repr__(self):
        rng = "%.1f ºC to %.1f ºC" % (self.func.minT, self.func.maxT)
        return "<%s thermocouple reference (%s)>" % (self.type, rng)

    @property
    def minT(self):
        return self.func.minT

    @property
    def maxT(self):
        return self.func.maxT

    def emfr(self, T, Tref=0.0, derivative=0, out_of_range="raise"):
        """
        Compute reference electromotive force for given thermocouple measurement 
        junction temperature and given reference junctions temperature.

        Parameters
        ----------
        T : array_like
            Temperature or array of temperatures (in ºC).
        Tref : float, optional
            Reference junctions' temperature (in ºC), defaults to 0.0.
            If derivative != 0, Tref is irrelevant.
        derivative : integer, optional
            Use this parameter to evaluate the functional derivative of
            the emf function at a given temperature.
            defaults to derivative=0 (no derivative).
        out_of_range : {'raise', 'nan', 'extrapolate'}, optional
            Determines behaviour for out of range temperatures: raise an
            exception, return NaNs, or extrapolate using the nearest
            polynomial. Note - do not trust the extrapolation!
        
        Returns
        -------
        emfr : array_like
            computed emfs (in mV)
            or, if derivative != 0,
            emf derivative (in mV / ºC**derivative)
        """
        f_T = self.func(T, derivative=derivative, out_of_range=out_of_range)
        if derivative != 0:
            return f_T
        f_ref = self.func(Tref, derivative=derivative, out_of_range=out_of_range)
        return f_T - f_ref

    def emfr_si(self, T, Tref=273.15, derivative=0, out_of_range="raise"):
        """
        This method is equivalent to emfr() but uses SI units - Kelvin and Volt.
        """
        return self.emfr(T - 273.15, Tref - 273.15, derivative, out_of_range) * 1.0e-3

    def t90r(self, emf, Tref=0.0, Tstart=None, Vtol=1.0e-6):
        """
        Inverse lookup: compute measurement junction temperature for a given
        measured voltage and given reference junctions temperature.
        
        Parameters
        ----------
        emf : float
            The measured voltage (in mV).
        Tref : float, optional
            The reference junctions temperature (in ºC), defaults to 0.0.
        Tstart : float, optional
            Suggested starting temperature (in ºC).
        Vtol : float, optional
            Tolerance of voltage in search, defaults to 1.0E-6 mV.
        
        Returns
        -------
        T : float
            Junction temperature (in ºC), such that:
              emf == func(T) - func(Tref)    (to within Vtol)
        """
        f_ref = self.func(Tref)
        T = self.func.inverse(emf + f_ref, Tstart=Tstart, Vtol=Vtol)
        return T

    def t90r_si(self, emf, Tref=273.15, Tstart=None, Vtol=1.0e-9):
        """
        This method is equivalent to t90r() but uses SI units - Kelvin and Volt.
        """
        return self.t90r(emf * 1.0e3, Tref - 273.15, Tstart, Vtol * 1.0e3)


# end of module
