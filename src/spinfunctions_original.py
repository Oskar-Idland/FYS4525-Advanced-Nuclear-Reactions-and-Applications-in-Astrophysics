import numpy as np
from scipy.interpolate import interp1d

"""
Code by Fabio Zeiser to calculate spin distribution
functions.
"""

class SpinFunctions:
    """ Calculates spin distributions, spin cuts (...)

    Args:
        Ex : double or dnarray
            Excitation energy
        J : double or dnarray
            Spin
        model : string
            Model to for the spincut
        pars : dict
            Additional parameters necessary for the spin cut model
    """

    def __init__(self, Ex, J=None, model=None, pars=None):
        self.Ex = np.atleast_1d(Ex)
        self.J = np.atleast_1d(J) if J is not None else None
        self.model = model
        self.pars = pars

    def get_sigma2(self):
        """ Get the square of the spin cut for a specified model """
        Ex = self.Ex
        model = self.model
        pars = self.pars

        assert model is not None and pars is not None

        def EB05(mass, NLDa, Eshift, Ex=Ex):
            Ex = np.atleast_1d(Ex)
            Eeff = Ex - Eshift
            Eeff[Eeff < 0] = 0
            sigma2 = (0.0146 * np.power(mass, 5.0 / 3.0)
                      * (1 + np.sqrt(1 + 4 * NLDa * Eeff))
                      / (2 * NLDa))
            return sigma2

        def EB09_CT(mass):
            sigma2 = np.power(0.98 * (mass**(0.29)), 2)
            return sigma2

        def EB09_emp(mass, Pa_prime, Ex=Ex):
            Ex = np.atleast_1d(Ex)
            Eeff = Ex - 0.5 * Pa_prime
            Eeff[Eeff < 0] = 0
            sigma2 = 0.391 * np.power(mass, 0.675) * np.power(Eeff, 0.312)
            return sigma2

        def Disc_and_EB05(mass, NLDa, Eshift, Sn, sigma2_disc):
            Ex_local = np.atleast_1d(Ex)
            sigma2_Sn = EB05(mass, NLDa, Eshift, Ex=np.atleast_1d(Sn))[0]
            sigma2_EB05_func = lambda E: EB05(mass, NLDa, Eshift, Ex=np.atleast_1d(E))[0]
            x = [sigma2_disc[0], Sn]
            y = [sigma2_disc[1], sigma2_EB05_func(Sn)]

            interp_sigma2 = interp1d(
                x, y,
                bounds_error=False,
                fill_value=(sigma2_disc[1], sigma2_Sn)
            )

            return np.where(Ex_local < Sn, interp_sigma2(Ex_local),
                            EB05(mass, NLDa, Eshift, Ex=Ex_local))

        if model == "EB05":
            pars_req = {"mass", "NLDa", "Eshift"}
            return call_model(EB05, pars, pars_req)
        elif model == "EB09_CT":
            pars_req = {"mass"}
            return call_model(EB09_CT, pars, pars_req)
        elif model == "EB09_emp":
            pars_req = {"mass", "Pa_prime"}
            return call_model(EB09_emp, pars, pars_req)
        elif model == "Disc_and_EB05":
            pars_req = {"mass", "NLDa", "Eshift", "Sn", "sigma2_disc"}
            return call_model(Disc_and_EB05, pars, pars_req)
        else:
            raise TypeError("\nError: Spincut model not supported; check spelling\n")

    def distribution(self):
        """ Get spin distribution assuming equal parity """
        J = self.J
        assert J is not None

        sigma2 = self.get_sigma2()
        sigma2 = sigma2[np.newaxis]

        spinDist = ((2. * J + 1.) / (2. * sigma2.T)
                    * np.exp(-np.power(J + 0.5, 2.) / (2. * sigma2.T)))
        return np.squeeze(spinDist)

def call_model(fun, pars, pars_req):
    """ Call `fun` and check if all required parameters are provided """
    if pars_req <= set(pars):
        pcall = {p: pars[p] for p in pars_req}
        return fun(**pcall)
    else:
        raise TypeError("Error: Need following arguments for this method: {0}".format(pars_req))
