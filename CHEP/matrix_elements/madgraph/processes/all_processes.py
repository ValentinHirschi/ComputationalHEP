from __future__ import division
from model.aloha_methods import *
from model.wavefunctions import *


class Matrix_2_epem_epem_no_z(object):

    def __init__(self):
        """define the object"""
        self.clean()

    def clean(self):
        self.jamp = []

    def get_external_masses(self, model):

        return ((model.ZERO, model.ZERO), (model.ZERO, model.ZERO))

    def smatrix(self, p, model):
        #
        #  MadGraph5_aMC@NLO v. 3.5.7, 2024-11-29
        #  By the MadGraph5_aMC@NLO Development Team
        #  Visit launchpad.net/madgraph5 and amcatnlo.web.cern.ch
        #
        # MadGraph5_aMC@NLO StandAlone Version
        #
        # Returns amplitude squared summed/avg over colors
        # and helicities
        # for the point in phase space P(0:3,NEXTERNAL)
        #
        # Process: e+ e- > e+ e- WEIGHTED<=4 / z @2
        #
        # Clean additional output
        #
        self.clean()
        #
        # CONSTANTS
        #
        nexternal = 4
        ndiags = 2
        ncomb = 16
        #
        # LOCAL VARIABLES
        #
        helicities = [
            [-1, 1, 1, -1],
            [-1, 1, 1, 1],
            [-1, 1, -1, -1],
            [-1, 1, -1, 1],
            [-1, -1, 1, -1],
            [-1, -1, 1, 1],
            [-1, -1, -1, -1],
            [-1, -1, -1, 1],
            [1, 1, 1, -1],
            [1, 1, 1, 1],
            [1, 1, -1, -1],
            [1, 1, -1, 1],
            [1, -1, 1, -1],
            [1, -1, 1, 1],
            [1, -1, -1, -1],
            [1, -1, -1, 1]]
        denominator = 4
        # ----------
        # BEGIN CODE
        # ----------
        self.amp2 = [0.] * ndiags
        self.helEvals = []
        ans = 0.
        for hel in helicities:
            t = self.matrix(p, hel, model)
            ans = ans + t
            self.helEvals.append([hel, t.real / denominator])
        ans = ans / denominator
        return ans.real

    def matrix(self, p, hel, model):
        #
        #  MadGraph5_aMC@NLO v. 3.5.7, 2024-11-29
        #  By the MadGraph5_aMC@NLO Development Team
        #  Visit launchpad.net/madgraph5 and amcatnlo.web.cern.ch
        #
        # Returns amplitude squared summed/avg over colors
        # for the point with external lines W(0:6,NEXTERNAL)
        #
        # Process: e+ e- > e+ e- WEIGHTED<=4 / z @2
        #
        #
        # Process parameters
        #
        ngraphs = 2
        nexternal = 4
        nwavefuncs = 5
        ncolor = 1
        ZERO = 0.
        #
        # Color matrix
        #
        denom = [1]
        cf = [[1]]
        #
        # Model parameters
        #

        GC_3 = model.GC_3
        # ----------
        # Begin code
        # ----------
        amp = [None] * ngraphs
        w = [None] * nwavefuncs
        w[0] = oxxxxx(p[0], ZERO, hel[0], -1)
        w[1] = ixxxxx(p[1], ZERO, hel[1], +1)
        w[2] = ixxxxx(p[2], ZERO, hel[2], -1)
        w[3] = oxxxxx(p[3], ZERO, hel[3], +1)
        w[4] = FFV1P0_3(w[1], w[0], GC_3, ZERO, ZERO)
        # Amplitude(s) for diagram number 1
        amp[0] = FFV1_0(w[2], w[3], w[4], GC_3)
        w[4] = FFV1P0_3(w[2], w[0], GC_3, ZERO, ZERO)
        # Amplitude(s) for diagram number 2
        amp[1] = FFV1_0(w[1], w[3], w[4], GC_3)

        jamp = [None] * ncolor

        jamp[0] = -amp[0]+amp[1]

        self.amp2[0] += abs(amp[0]*amp[0].conjugate())
        self.amp2[1] += abs(amp[1]*amp[1].conjugate())
        matrix = 0.
        for i in range(ncolor):
            ztemp = 0
            for j in range(ncolor):
                ztemp = ztemp + cf[i][j]*jamp[j]
            matrix = matrix + ztemp * jamp[i].conjugate()/denom[i]
        self.jamp.append(jamp)

        return matrix


class Matrix_3_epem_ddxg_no_z(object):

    def __init__(self):
        """define the object"""
        self.clean()

    def clean(self):
        self.jamp = []

    def get_external_masses(self, model):

        return ((model.ZERO, model.ZERO), (model.ZERO, model.ZERO, model.ZERO))

    def smatrix(self, p, model):
        #
        #  MadGraph5_aMC@NLO v. 3.5.7, 2024-11-29
        #  By the MadGraph5_aMC@NLO Development Team
        #  Visit launchpad.net/madgraph5 and amcatnlo.web.cern.ch
        #
        # MadGraph5_aMC@NLO StandAlone Version
        #
        # Returns amplitude squared summed/avg over colors
        # and helicities
        # for the point in phase space P(0:3,NEXTERNAL)
        #
        # Process: e+ e- > d d~ g WEIGHTED<=5 / z @3
        #
        # Clean additional output
        #
        self.clean()
        #
        # CONSTANTS
        #
        nexternal = 5
        ndiags = 2
        ncomb = 32
        #
        # LOCAL VARIABLES
        #
        helicities = [
            [-1, 1, -1, 1, -1],
            [-1, 1, -1, 1, 1],
            [-1, 1, -1, -1, -1],
            [-1, 1, -1, -1, 1],
            [-1, 1, 1, 1, -1],
            [-1, 1, 1, 1, 1],
            [-1, 1, 1, -1, -1],
            [-1, 1, 1, -1, 1],
            [-1, -1, -1, 1, -1],
            [-1, -1, -1, 1, 1],
            [-1, -1, -1, -1, -1],
            [-1, -1, -1, -1, 1],
            [-1, -1, 1, 1, -1],
            [-1, -1, 1, 1, 1],
            [-1, -1, 1, -1, -1],
            [-1, -1, 1, -1, 1],
            [1, 1, -1, 1, -1],
            [1, 1, -1, 1, 1],
            [1, 1, -1, -1, -1],
            [1, 1, -1, -1, 1],
            [1, 1, 1, 1, -1],
            [1, 1, 1, 1, 1],
            [1, 1, 1, -1, -1],
            [1, 1, 1, -1, 1],
            [1, -1, -1, 1, -1],
            [1, -1, -1, 1, 1],
            [1, -1, -1, -1, -1],
            [1, -1, -1, -1, 1],
            [1, -1, 1, 1, -1],
            [1, -1, 1, 1, 1],
            [1, -1, 1, -1, -1],
            [1, -1, 1, -1, 1]]
        denominator = 4
        # ----------
        # BEGIN CODE
        # ----------
        self.amp2 = [0.] * ndiags
        self.helEvals = []
        ans = 0.
        for hel in helicities:
            t = self.matrix(p, hel, model)
            ans = ans + t
            self.helEvals.append([hel, t.real / denominator])
        ans = ans / denominator
        return ans.real

    def matrix(self, p, hel, model):
        #
        #  MadGraph5_aMC@NLO v. 3.5.7, 2024-11-29
        #  By the MadGraph5_aMC@NLO Development Team
        #  Visit launchpad.net/madgraph5 and amcatnlo.web.cern.ch
        #
        # Returns amplitude squared summed/avg over colors
        # for the point with external lines W(0:6,NEXTERNAL)
        #
        # Process: e+ e- > d d~ g WEIGHTED<=5 / z @3
        #
        #
        # Process parameters
        #
        ngraphs = 2
        nexternal = 5
        nwavefuncs = 6
        ncolor = 1
        ZERO = 0.
        #
        # Color matrix
        #
        denom = [1]
        cf = [[4]]
        #
        # Model parameters
        #

        GC_3 = model.GC_3
        GC_11 = model.GC_11
        GC_1 = model.GC_1
        # ----------
        # Begin code
        # ----------
        amp = [None] * ngraphs
        w = [None] * nwavefuncs
        w[0] = oxxxxx(p[0], ZERO, hel[0], -1)
        w[1] = ixxxxx(p[1], ZERO, hel[1], +1)
        w[2] = oxxxxx(p[2], ZERO, hel[2], +1)
        w[3] = ixxxxx(p[3], ZERO, hel[3], -1)
        w[4] = vxxxxx(p[4], ZERO, hel[4], +1)
        w[5] = FFV1P0_3(w[1], w[0], GC_3, ZERO, ZERO)
        w[1] = FFV1_1(w[2], w[4], GC_11, ZERO, ZERO)
        # Amplitude(s) for diagram number 1
        amp[0] = FFV1_0(w[3], w[1], w[5], GC_1)
        w[1] = FFV1_2(w[3], w[4], GC_11, ZERO, ZERO)
        # Amplitude(s) for diagram number 2
        amp[1] = FFV1_0(w[1], w[2], w[5], GC_1)

        jamp = [None] * ncolor

        jamp[0] = +amp[0]+amp[1]

        self.amp2[0] += abs(amp[0]*amp[0].conjugate())
        self.amp2[1] += abs(amp[1]*amp[1].conjugate())
        matrix = 0.
        for i in range(ncolor):
            ztemp = 0
            for j in range(ncolor):
                ztemp = ztemp + cf[i][j]*jamp[j]
            matrix = matrix + ztemp * jamp[i].conjugate()/denom[i]
        self.jamp.append(jamp)

        return matrix


class Matrix_1_epem_mupmum_no_z(object):

    def __init__(self):
        """define the object"""
        self.clean()

    def clean(self):
        self.jamp = []

    def get_external_masses(self, model):

        return ((model.ZERO, model.ZERO), (model.ZERO, model.ZERO))

    def smatrix(self, p, model):
        #
        #  MadGraph5_aMC@NLO v. 3.5.7, 2024-11-29
        #  By the MadGraph5_aMC@NLO Development Team
        #  Visit launchpad.net/madgraph5 and amcatnlo.web.cern.ch
        #
        # MadGraph5_aMC@NLO StandAlone Version
        #
        # Returns amplitude squared summed/avg over colors
        # and helicities
        # for the point in phase space P(0:3,NEXTERNAL)
        #
        # Process: e+ e- > mu+ mu- WEIGHTED<=4 / z @1
        #
        # Clean additional output
        #
        self.clean()
        #
        # CONSTANTS
        #
        nexternal = 4
        ndiags = 1
        ncomb = 16
        #
        # LOCAL VARIABLES
        #
        helicities = [
            [-1, 1, 1, -1],
            [-1, 1, 1, 1],
            [-1, 1, -1, -1],
            [-1, 1, -1, 1],
            [-1, -1, 1, -1],
            [-1, -1, 1, 1],
            [-1, -1, -1, -1],
            [-1, -1, -1, 1],
            [1, 1, 1, -1],
            [1, 1, 1, 1],
            [1, 1, -1, -1],
            [1, 1, -1, 1],
            [1, -1, 1, -1],
            [1, -1, 1, 1],
            [1, -1, -1, -1],
            [1, -1, -1, 1]]
        denominator = 4
        # ----------
        # BEGIN CODE
        # ----------
        self.amp2 = [0.] * ndiags
        self.helEvals = []
        ans = 0.
        for hel in helicities:
            t = self.matrix(p, hel, model)
            ans = ans + t
            self.helEvals.append([hel, t.real / denominator])
        ans = ans / denominator
        return ans.real

    def matrix(self, p, hel, model):
        #
        #  MadGraph5_aMC@NLO v. 3.5.7, 2024-11-29
        #  By the MadGraph5_aMC@NLO Development Team
        #  Visit launchpad.net/madgraph5 and amcatnlo.web.cern.ch
        #
        # Returns amplitude squared summed/avg over colors
        # for the point with external lines W(0:6,NEXTERNAL)
        #
        # Process: e+ e- > mu+ mu- WEIGHTED<=4 / z @1
        #
        #
        # Process parameters
        #
        ngraphs = 1
        nexternal = 4
        nwavefuncs = 5
        ncolor = 1
        ZERO = 0.
        #
        # Color matrix
        #
        denom = [1]
        cf = [[1]]
        #
        # Model parameters
        #

        GC_3 = model.GC_3
        # ----------
        # Begin code
        # ----------
        amp = [None] * ngraphs
        w = [None] * nwavefuncs
        w[0] = oxxxxx(p[0], ZERO, hel[0], -1)
        w[1] = ixxxxx(p[1], ZERO, hel[1], +1)
        w[2] = ixxxxx(p[2], ZERO, hel[2], -1)
        w[3] = oxxxxx(p[3], ZERO, hel[3], +1)
        w[4] = FFV1P0_3(w[1], w[0], GC_3, ZERO, ZERO)
        # Amplitude(s) for diagram number 1
        amp[0] = FFV1_0(w[2], w[3], w[4], GC_3)

        jamp = [None] * ncolor

        jamp[0] = -amp[0]

        self.amp2[0] += abs(amp[0]*amp[0].conjugate())
        matrix = 0.
        for i in range(ncolor):
            ztemp = 0
            for j in range(ncolor):
                ztemp = ztemp + cf[i][j]*jamp[j]
            matrix = matrix + ztemp * jamp[i].conjugate()/denom[i]
        self.jamp.append(jamp)

        return matrix


class Matrix_1_uux_wpwm_no_za(object):

    def __init__(self):
        """define the object"""
        self.clean()

    def clean(self):
        self.jamp = []

    def get_external_masses(self, model):

        return ((model.ZERO, model.ZERO), (model.mdl_MW, model.mdl_MW))

    def smatrix(self, p, model):
        #
        #  MadGraph5_aMC@NLO v. 3.5.7, 2024-11-29
        #  By the MadGraph5_aMC@NLO Development Team
        #  Visit launchpad.net/madgraph5 and amcatnlo.web.cern.ch
        #
        # MadGraph5_aMC@NLO StandAlone Version
        #
        # Returns amplitude squared summed/avg over colors
        # and helicities
        # for the point in phase space P(0:3,NEXTERNAL)
        #
        # Process: u u~ > w+ w- WEIGHTED<=4 / z a @1
        # Process: c c~ > w+ w- WEIGHTED<=4 / z a @1
        #
        # Clean additional output
        #
        self.clean()
        #
        # CONSTANTS
        #
        nexternal = 4
        ndiags = 1
        ncomb = 36
        #
        # LOCAL VARIABLES
        #
        helicities = [
            [1, -1, -1, 1],
            [1, -1, -1, 0],
            [1, -1, -1, -1],
            [1, -1, 0, 1],
            [1, -1, 0, 0],
            [1, -1, 0, -1],
            [1, -1, 1, 1],
            [1, -1, 1, 0],
            [1, -1, 1, -1],
            [1, 1, -1, 1],
            [1, 1, -1, 0],
            [1, 1, -1, -1],
            [1, 1, 0, 1],
            [1, 1, 0, 0],
            [1, 1, 0, -1],
            [1, 1, 1, 1],
            [1, 1, 1, 0],
            [1, 1, 1, -1],
            [-1, -1, -1, 1],
            [-1, -1, -1, 0],
            [-1, -1, -1, -1],
            [-1, -1, 0, 1],
            [-1, -1, 0, 0],
            [-1, -1, 0, -1],
            [-1, -1, 1, 1],
            [-1, -1, 1, 0],
            [-1, -1, 1, -1],
            [-1, 1, -1, 1],
            [-1, 1, -1, 0],
            [-1, 1, -1, -1],
            [-1, 1, 0, 1],
            [-1, 1, 0, 0],
            [-1, 1, 0, -1],
            [-1, 1, 1, 1],
            [-1, 1, 1, 0],
            [-1, 1, 1, -1]]
        denominator = 36
        # ----------
        # BEGIN CODE
        # ----------
        self.amp2 = [0.] * ndiags
        self.helEvals = []
        ans = 0.
        for hel in helicities:
            t = self.matrix(p, hel, model)
            ans = ans + t
            self.helEvals.append([hel, t.real / denominator])
        ans = ans / denominator
        return ans.real

    def matrix(self, p, hel, model):
        #
        #  MadGraph5_aMC@NLO v. 3.5.7, 2024-11-29
        #  By the MadGraph5_aMC@NLO Development Team
        #  Visit launchpad.net/madgraph5 and amcatnlo.web.cern.ch
        #
        # Returns amplitude squared summed/avg over colors
        # for the point with external lines W(0:6,NEXTERNAL)
        #
        # Process: u u~ > w+ w- WEIGHTED<=4 / z a @1
        # Process: c c~ > w+ w- WEIGHTED<=4 / z a @1
        #
        #
        # Process parameters
        #
        ngraphs = 1
        nexternal = 4
        nwavefuncs = 5
        ncolor = 1
        ZERO = 0.
        #
        # Color matrix
        #
        denom = [1]
        cf = [[3]]
        #
        # Model parameters
        #
        mdl_WW = model.mdl_WW
        mdl_MW = model.mdl_MW
        GC_100 = model.GC_100
        # ----------
        # Begin code
        # ----------
        amp = [None] * ngraphs
        w = [None] * nwavefuncs
        w[0] = ixxxxx(p[0], ZERO, hel[0], +1)
        w[1] = oxxxxx(p[1], ZERO, hel[1], -1)
        w[2] = vxxxxx(p[2], mdl_MW, hel[2], +1)
        w[3] = vxxxxx(p[3], mdl_MW, hel[3], +1)
        w[4] = FFV2_2(w[0], w[2], GC_100, ZERO, ZERO)
        # Amplitude(s) for diagram number 1
        amp[0] = FFV2_0(w[4], w[1], w[3], GC_100)

        jamp = [None] * ncolor

        jamp[0] = +amp[0]

        self.amp2[0] += abs(amp[0]*amp[0].conjugate())
        matrix = 0.
        for i in range(ncolor):
            ztemp = 0
            for j in range(ncolor):
                ztemp = ztemp + cf[i][j]*jamp[j]
            matrix = matrix + ztemp * jamp[i].conjugate()/denom[i]
        self.jamp.append(jamp)

        return matrix


class Matrix_1_ddx_wpwm_no_za(object):

    def __init__(self):
        """define the object"""
        self.clean()

    def clean(self):
        self.jamp = []

    def get_external_masses(self, model):

        return ((model.ZERO, model.ZERO), (model.mdl_MW, model.mdl_MW))

    def smatrix(self, p, model):
        #
        #  MadGraph5_aMC@NLO v. 3.5.7, 2024-11-29
        #  By the MadGraph5_aMC@NLO Development Team
        #  Visit launchpad.net/madgraph5 and amcatnlo.web.cern.ch
        #
        # MadGraph5_aMC@NLO StandAlone Version
        #
        # Returns amplitude squared summed/avg over colors
        # and helicities
        # for the point in phase space P(0:3,NEXTERNAL)
        #
        # Process: d d~ > w+ w- WEIGHTED<=4 / z a @1
        # Process: s s~ > w+ w- WEIGHTED<=4 / z a @1
        #
        # Clean additional output
        #
        self.clean()
        #
        # CONSTANTS
        #
        nexternal = 4
        ndiags = 1
        ncomb = 36
        #
        # LOCAL VARIABLES
        #
        helicities = [
            [1, -1, -1, 1],
            [1, -1, -1, 0],
            [1, -1, -1, -1],
            [1, -1, 0, 1],
            [1, -1, 0, 0],
            [1, -1, 0, -1],
            [1, -1, 1, 1],
            [1, -1, 1, 0],
            [1, -1, 1, -1],
            [1, 1, -1, 1],
            [1, 1, -1, 0],
            [1, 1, -1, -1],
            [1, 1, 0, 1],
            [1, 1, 0, 0],
            [1, 1, 0, -1],
            [1, 1, 1, 1],
            [1, 1, 1, 0],
            [1, 1, 1, -1],
            [-1, -1, -1, 1],
            [-1, -1, -1, 0],
            [-1, -1, -1, -1],
            [-1, -1, 0, 1],
            [-1, -1, 0, 0],
            [-1, -1, 0, -1],
            [-1, -1, 1, 1],
            [-1, -1, 1, 0],
            [-1, -1, 1, -1],
            [-1, 1, -1, 1],
            [-1, 1, -1, 0],
            [-1, 1, -1, -1],
            [-1, 1, 0, 1],
            [-1, 1, 0, 0],
            [-1, 1, 0, -1],
            [-1, 1, 1, 1],
            [-1, 1, 1, 0],
            [-1, 1, 1, -1]]
        denominator = 36
        # ----------
        # BEGIN CODE
        # ----------
        self.amp2 = [0.] * ndiags
        self.helEvals = []
        ans = 0.
        for hel in helicities:
            t = self.matrix(p, hel, model)
            ans = ans + t
            self.helEvals.append([hel, t.real / denominator])
        ans = ans / denominator
        return ans.real

    def matrix(self, p, hel, model):
        #
        #  MadGraph5_aMC@NLO v. 3.5.7, 2024-11-29
        #  By the MadGraph5_aMC@NLO Development Team
        #  Visit launchpad.net/madgraph5 and amcatnlo.web.cern.ch
        #
        # Returns amplitude squared summed/avg over colors
        # for the point with external lines W(0:6,NEXTERNAL)
        #
        # Process: d d~ > w+ w- WEIGHTED<=4 / z a @1
        # Process: s s~ > w+ w- WEIGHTED<=4 / z a @1
        #
        #
        # Process parameters
        #
        ngraphs = 1
        nexternal = 4
        nwavefuncs = 5
        ncolor = 1
        ZERO = 0.
        #
        # Color matrix
        #
        denom = [1]
        cf = [[3]]
        #
        # Model parameters
        #
        mdl_WW = model.mdl_WW
        mdl_MW = model.mdl_MW
        GC_100 = model.GC_100
        # ----------
        # Begin code
        # ----------
        amp = [None] * ngraphs
        w = [None] * nwavefuncs
        w[0] = ixxxxx(p[0], ZERO, hel[0], +1)
        w[1] = oxxxxx(p[1], ZERO, hel[1], -1)
        w[2] = vxxxxx(p[2], mdl_MW, hel[2], +1)
        w[3] = vxxxxx(p[3], mdl_MW, hel[3], +1)
        w[4] = FFV2_2(w[0], w[3], GC_100, ZERO, ZERO)
        # Amplitude(s) for diagram number 1
        amp[0] = FFV2_0(w[4], w[1], w[2], GC_100)

        jamp = [None] * ncolor

        jamp[0] = +amp[0]

        self.amp2[0] += abs(amp[0]*amp[0].conjugate())
        matrix = 0.
        for i in range(ncolor):
            ztemp = 0
            for j in range(ncolor):
                ztemp = ztemp + cf[i][j]*jamp[j]
            matrix = matrix + ztemp * jamp[i].conjugate()/denom[i]
        self.jamp.append(jamp)

        return matrix
