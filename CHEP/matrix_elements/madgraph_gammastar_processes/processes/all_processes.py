from __future__ import division
from model.aloha_methods import *
from model.wavefunctions import *
class Matrix_2_z_ddxg(object):

    def __init__(self):
        """define the object"""
        self.clean()

    def clean(self):
        self.jamp = []

    def get_external_masses(self, model):

        return ( (model.mdl_MZ), (model.ZERO, model.ZERO, model.ZERO) )

    def smatrix(self,p, model):
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
        # Process: z > d d~ g @2
        #  
        # Clean additional output
        #
        self.clean()
        #  
        # CONSTANTS
        #  
        nexternal = 4
        ndiags = 2
        ncomb = 24
        #  
        # LOCAL VARIABLES 
        #  
        helicities = [ \
        [-1,-1,1,-1],
        [-1,-1,1,1],
        [-1,-1,-1,-1],
        [-1,-1,-1,1],
        [-1,1,1,-1],
        [-1,1,1,1],
        [-1,1,-1,-1],
        [-1,1,-1,1],
        [0,-1,1,-1],
        [0,-1,1,1],
        [0,-1,-1,-1],
        [0,-1,-1,1],
        [0,1,1,-1],
        [0,1,1,1],
        [0,1,-1,-1],
        [0,1,-1,1],
        [1,-1,1,-1],
        [1,-1,1,1],
        [1,-1,-1,-1],
        [1,-1,-1,1],
        [1,1,1,-1],
        [1,1,1,1],
        [1,1,-1,-1],
        [1,1,-1,1]]
        denominator = 3
        # ----------
        # BEGIN CODE
        # ----------
        self.amp2 = [0.] * ndiags
        self.helEvals = []
        ans = 0.
        for hel in helicities:
            t = self.matrix(p, hel, model)
            ans = ans + t
            self.helEvals.append([hel, t.real / denominator ])
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
        # Process: z > d d~ g @2
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
        denom = [1];
        cf = [[4]];
        #
        # Model parameters
        #
        mdl_MZ = model.mdl_MZ
        mdl_WZ = model.mdl_WZ
        GC_1 = model.GC_1
        GC_11 = model.GC_11
        # ----------
        # Begin code
        # ----------
        amp = [None] * ngraphs
        w = [None] * nwavefuncs
        w[0] = vxxxxx(p[0],mdl_MZ,hel[0],-1)
        w[1] = oxxxxx(p[1],ZERO,hel[1],+1)
        w[2] = ixxxxx(p[2],ZERO,hel[2],-1)
        w[3] = vxxxxx(p[3],ZERO,hel[3],+1)
        w[4]= FFV1_1(w[1],w[3],GC_11,ZERO,ZERO)
        # Amplitude(s) for diagram number 1
        amp[0]= FFV1_0(w[2],w[4],w[0],GC_1)
        w[4]= FFV1_2(w[2],w[3],GC_11,ZERO,ZERO)
        # Amplitude(s) for diagram number 2
        amp[1]= FFV1_0(w[4],w[1],w[0],GC_1)

        jamp = [None] * ncolor

        jamp[0] = -amp[0]-amp[1]

        self.amp2[0]+=abs(amp[0]*amp[0].conjugate())
        self.amp2[1]+=abs(amp[1]*amp[1].conjugate())
        matrix = 0.
        for i in range(ncolor):
            ztemp = 0
            for j in range(ncolor):
                ztemp = ztemp + cf[i][j]*jamp[j]
            matrix = matrix + ztemp * jamp[i].conjugate()/denom[i]   
        self.jamp.append(jamp)

        return matrix

class Matrix_1_z_ddx(object):

    def __init__(self):
        """define the object"""
        self.clean()

    def clean(self):
        self.jamp = []

    def get_external_masses(self, model):

        return ( (model.mdl_MZ), (model.ZERO, model.ZERO) )

    def smatrix(self,p, model):
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
        # Process: z > d d~ @1
        #  
        # Clean additional output
        #
        self.clean()
        #  
        # CONSTANTS
        #  
        nexternal = 3
        ndiags = 1
        ncomb = 12
        #  
        # LOCAL VARIABLES 
        #  
        helicities = [ \
        [-1,-1,1],
        [-1,-1,-1],
        [-1,1,1],
        [-1,1,-1],
        [0,-1,1],
        [0,-1,-1],
        [0,1,1],
        [0,1,-1],
        [1,-1,1],
        [1,-1,-1],
        [1,1,1],
        [1,1,-1]]
        denominator = 3
        # ----------
        # BEGIN CODE
        # ----------
        self.amp2 = [0.] * ndiags
        self.helEvals = []
        ans = 0.
        for hel in helicities:
            t = self.matrix(p, hel, model)
            ans = ans + t
            self.helEvals.append([hel, t.real / denominator ])
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
        # Process: z > d d~ @1
        #  
        #  
        # Process parameters
        #  
        ngraphs = 1
        nexternal = 3
        nwavefuncs = 3
        ncolor = 1
        ZERO = 0.
        #  
        # Color matrix
        #  
        denom = [1];
        cf = [[3]];
        #
        # Model parameters
        #
        mdl_MZ = model.mdl_MZ
        mdl_WZ = model.mdl_WZ
        GC_1 = model.GC_1
        # ----------
        # Begin code
        # ----------
        amp = [None] * ngraphs
        w = [None] * nwavefuncs
        w[0] = vxxxxx(p[0],mdl_MZ,hel[0],-1)
        w[1] = oxxxxx(p[1],ZERO,hel[1],+1)
        w[2] = ixxxxx(p[2],ZERO,hel[2],-1)
        # Amplitude(s) for diagram number 1
        amp[0]= FFV1_0(w[2],w[1],w[0],GC_1)

        jamp = [None] * ncolor

        jamp[0] = -amp[0]

        self.amp2[0]+=abs(amp[0]*amp[0].conjugate())
        matrix = 0.
        for i in range(ncolor):
            ztemp = 0
            for j in range(ncolor):
                ztemp = ztemp + cf[i][j]*jamp[j]
            matrix = matrix + ztemp * jamp[i].conjugate()/denom[i]   
        self.jamp.append(jamp)

        return matrix

